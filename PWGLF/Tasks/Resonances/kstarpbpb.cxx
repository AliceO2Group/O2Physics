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

/// \file kstarpbpb.cxx
/// \brief Code for K*(892)^0 resonance flow and spin alignment analysis
/// \author sourav.kundu@cern.ch , sarjeeta.gami@cern.ch
///

#include "PWGLF/DataModel/EPCalibrationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
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
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <THn.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct Kstarpbpb {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"nolaterthan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::ccdb::CcdbApi ccdbApi;
  // Service<o2::framework::O2DatabasePDG> pdg;
  struct RCTCut : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};

    RCTFlagsChecker rctChecker;
  };

  RCTCut rctCut;

  // CCDB options
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  // Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  // track
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<bool> usepolar{"usepolar", true, "flag to fill type of SA"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<bool> isTOFOnly{"isTOFOnly", false, "use TOF only PID"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {180, 0.6, 1.5}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configrapAxis{"configrapAxis", {VARIABLE_WIDTH, -0.8, -0.4, 0.4, 0.8}, "Rapidity"};
  Configurable<float> confFakeKaonCut{"confFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {400, -16, 16}, "V2"};
  Configurable<bool> isNoTOF{"isNoTOF", true, "isNoTOF"};
  Configurable<bool> pdgcheck{"pdgcheck", true, "pdgcheck"};
  Configurable<int> strategyPID{"strategyPID", 2, "PID strategy"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalQAplots1{"additionalQAplots1", true, "Additional QA plots"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * o2::constants::math::PI / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * o2::constants::math::PI / 6.0, "Maximum of rotation"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<bool> fillSA{"fillSA", true, "same event SA"};
  Configurable<bool> doClosureTest{"doClosureTest", false, "Fill dedicated MC closure histograms for spin alignment"};
  Configurable<bool> fillLikeSignSA{"fillLikeSignSA", false, "fill the same-event like-sign SA sparse (needs fillSA)"};
  Configurable<bool> useWeight{"useWeight", false, "use EP dep effi weight"};
  Configurable<bool> useSP{"useSP", false, "use SP"};
  // spin-alignment quantization axis
  Configurable<int> cfgSAFrame{"cfgSAFrame", 0, "SA quantization axis: 0 = event plane (FT0C), 1 = production plane, 2 = random event plane"};
  Configurable<int> cfgRndSeed{"cfgRndSeed", 0, "Seed for random event plane (0 = unique seed per job)"};
  Configurable<float> cfgEPNormalHarmonic{"cfgEPNormalHarmonic", 1.0f, "EP normal = (sin(n*Psi), -cos(n*Psi), 0); n = 1 is the normal to the Psi2 plane (default), n = 2 is the old behaviour (axis at 2*Psi2 - pi/2, not the normal; the R2 resolution correction does not apply)"};
  // phi(1020) -> K+K- reflection veto for K* (off by default: the task then behaves exactly as before)
  Configurable<bool> cfgPhiVeto{"cfgPhiVeto", false, "Reject K pi pairs whose KK mass (pion candidate given the kaon mass) is near the phi(1020) mass"};
  Configurable<float> cfgPhiVetoWindow{"cfgPhiVetoWindow", 0.010f, "Half-width of the phi(1020) veto window (GeV/c^2)"};

  // event loss and signal loss for K* and phi(1020) (processEvtLossSigLossMC), adapted from processEvtLossSigLossMC1 in phianalysisrun3pbpb.cxx
  struct : ConfigurableGroup {
    std::string prefix = "evtSigLoss";
    Configurable<bool> cutVzGen{"cutVzGen", true, "Cut |vz| < cfgCutVertex on the generated collision"};
    Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Require INEL > 0 for the generated collision"};
    Configurable<bool> isApplyTVX{"isApplyTVX", false, "Require TVX (FT0A and FT0C MC multiplicity > 0) for the generated collision"};
    Configurable<bool> requireSingleReco{"requireSingleReco", true, "Accepted MC event: exactly one reconstructed collision passing the selection, as in processMC / processMCPhi (false: at least one, as in phianalysisrun3pbpb.cxx)"};
    ConfigurableAxis axisImpactPar{"axisImpactPar", {200, 0.0, 20.0}, "Impact parameter (fm)"};
    ConfigurableAxis axisMultMC{"axisMultMC", {500, 0.0, 2000.0}, "Generated charged particles in |eta| < 0.5"};
  } evtSigLoss;

  // phi(1020) spin alignment: event selection is the K* one; track and PID selections identical to PWGLF/Tasks/Resonances/phianalysisrun3pbpb.cxx
  // (JSON keys follow the O2 linter naming rules)
  struct : ConfigurableGroup {
    std::string prefix = "phiSA";
    Configurable<bool> removeFakeTrack{"removeFakeTrack", true, "Remove fake track from momentum difference"};
    Configurable<float> confFakeKaonCut{"confFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
    // track selection (phianalysisrun3pbpb.cxx)
    Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "Global track (with DCA) + PV contributor + ITS clusters"};
    Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "Global track w/o DCA + PV contributor + manual DCAxy/DCAz + ITS clusters"};
    Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for kaon tracks (manual DCA cut)"};
    Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for kaon tracks (manual DCA cut)"};
    Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
    // PID (phianalysisrun3pbpb.cxx)
    Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
    Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the combined TPC-TOF Nsigma cut"};
    Configurable<bool> isNoTOF{"isNoTOF", false, "TPC-only PID for all tracks"};
    Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
    Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
    Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
    Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
    Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
    ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
    ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, -1.0, 1.}, "cos(#vartheta)"};
    ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
    ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {8, 0, 0.8}, "Rapidity"};
    ConfigurableAxis configThnAxisSA{"configThnAxisSA", {200, -1, 1}, "SA"};
  } phiSA;

  // MC closure test of the cos(theta*) chain (processMCClosure, processMCClosureME, processMCClosurePhi, processMCClosureMEPhi)
  // Step 1: SE unlike-sign minus normalised ME (same cuts as data) -> raw yield vs cos(theta*)
  // Step 2: raw yield compared with the true reconstructed yield (hSparseRecTrue)
  // Step 3: raw yield / (Acc x Eff)(cos theta*), with Acc x Eff = hSparseRecTrueEff / hSparseGenEff
  // Step 4: corrected yield and rho00 compared with the generated input (hSparseGen)
  struct : ConfigurableGroup {
    std::string prefix = "closure";
    Configurable<bool> injectPolarization{"injectPolarization", false, "Reweight the true K*/phi with W(cos theta*) so that the sample has rho00 = rho00Inj"};
    Configurable<float> rho00Inj{"rho00Inj", 1.0f / 3.0f, "Injected rho00, in [0,1] (1/3 = unpolarized)"};
    Configurable<bool> useTrueEP{"useTrueEP", true, "Use the MC true event plane angle as Psi (false: Psi = 0, fixed axis as in processMC)"};
    Configurable<bool> splitSample{"splitSample", true, "Even MC events build Acc x Eff, odd MC events are the test sample (false: same events for both)"};
    Configurable<bool> fillRotBkg{"fillRotBkg", false, "Fill the rotated-pair background (weight 1), cross-check of the ME background"};
    Configurable<bool> genAcceptanceCut{"genAcceptanceCut", false, "Daughter pT/eta cut on the generated K*/phi (false: denominator is the full |y| < cut sample, so the correction is the full Acc x Eff in cos theta*)"};
    Configurable<bool> requirePrimaryTracks{"requirePrimaryTracks", false, "SE, ME and rotated pairs only from tracks matched to physical primaries (as in phianalysisrun3pbpb / kstar892LightIon); false: all tracks, as in data"};
    Configurable<bool> meRequireSingleReco{"meRequireSingleReco", true, "ME: both events must come from an MC collision with exactly one reco collision, as required in the SE (processMCClosure / processMCClosurePhi)"};
    ConfigurableAxis axisCosThetaStar{"axisCosThetaStar", {10, -1.0, 1.0}, "cos(#vartheta*)"};
  } closure;
  Configurable<float> cfgMinTrackPt{"cfgMinTrackPt", 0.15f,
                                    "Minimum track pT"};

  Configurable<float> cfgMaxTrackPt{"cfgMaxTrackPt", 10.0f,
                                    "Maximum track pT"};
  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<bool> additionalEvSel1{"additionalEvSel1", true, "Additional evsel1"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<bool> additionalEvSel4{"additionalEvSel4", true, "Additional evsel4"};
  Configurable<std::string> confWeightPath{"confWeightPath", "Users/s/skundu/My/Object/fitweight", "Path to gain calibration"};
  ConfigurableAxis axisPtKaonWeight{"axisPtKaonWeight", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcacutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  using CollisionMCTrueTable = aod::McCollisions;
  using McCollisionMults = soa::Join<aod::McCollisions, aod::MultMCExtras>;
  using TrackMCTrueTable = aod::McParticles;
  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;
  // reconstructed MC collisions for the closure: adds EPCalibrationTables (psiFT0C etc.) on top of CollisionMCRecTableCentFT0C,
  // so the closure can optionally use the reconstructed (finite-resolution) Psi_FT0C instead of the MC true event plane
  using CollisionMCRecTableClosure = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::EPCalibrationTables>>;
  // reconstructed MC collisions for the closure mixed event (posZ and centrality filters as in data)
  using CollisionMCRecMixTable = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs, aod::EvSels>>;
  // same as CollisionMCRecMixTable, plus EPCalibrationTables (psiFT0C etc.) for the reconstructed-event-plane closure
  using CollisionMCRecMixTableRecoEP = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs, aod::EvSels, aod::EPCalibrationTables>>;
  using BinningTypeMCClosure = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  // same as BinningTypeMCClosure, plus Psi_FT0C (as in data's BinningTypeVertexContributor), for the reconstructed-EP closure ME
  using BinningTypeMCClosureRecoEP = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;

  SliceCache cache;
  // Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  // Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  enum SAFrame : int {
    kEventPlane = 0,
    kProductionPlane = 1,
    kRandomEventPlane = 2
  };
  TRandom3 rndGen;

  void init(o2::framework::InitContext&)
  {
    rctCut.rctChecker.init(
      rctCut.cfgEvtRCTFlagCheckerLabel,
      rctCut.cfgEvtRCTFlagCheckerZDCCheck,
      rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    if (cfgSAFrame.value < kEventPlane || cfgSAFrame.value > kRandomEventPlane) {
      LOGF(fatal, "cfgSAFrame = %d not supported (0 = event plane, 1 = production plane, 2 = random event plane)", cfgSAFrame.value);
    }
    if (cfgSAFrame.value != kProductionPlane && cfgEPNormalHarmonic.value != 1.0f) {
      LOGF(warning, "cfgEPNormalHarmonic = %.1f: the SA axis is not the normal to the Psi2 plane, so the R2 = <cos 2(Psi_obs - Psi_true)> resolution correction does not close", cfgEPNormalHarmonic.value);
    }
    rndGen.SetSeed(cfgRndSeed);

    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {6000, -30, 30, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    AxisSpec vzAxis = {400, -20.0, 20.0, "Z_{vtx} (cm)"};
    AxisSpec qvecAxis = {200, 0, 20, "q"};
    if (doprocessSE) {
      histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{10, 0, 10.0}});
      if (!fillSA) {
        histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      }
      if (fillRotation && !fillSA) {
        histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});
        histos.add("hSparseV2SASameEventRotational_V2", "hSparseV2SASameEventRotational_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      }
      if (fillSA) {
        histos.add("hSparseSAvsrapsameunlike", "hSparseSAvsrapsameunlike", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
        if (fillLikeSignSA) {
          histos.add("hSparseSAvsrapsamelike", "hSparseSAvsrapsamelike", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
        }
        if (fillRotation) {
          histos.add("hSparseSAvsraprot", "hSparseSAvsraprot", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
        }
        if (cfgSAFrame.value == kRandomEventPlane) {
          histos.add("hPsiRandom", "Random event plane angle", kTH2F, {centAxis, phiAxis});
        }
      }
    }

    if (doprocessMixedEvent) {
      if (!fillSA) {
        histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      }
      if (fillSA) {
        histos.add("hSparseSAvsrapmix", "hSparseSAvsrapmix", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
      }
    }

    if (doprocessMC) {
      histos.add("hSparseV2SAGen_V2", "hSparseV2SAGen_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      histos.add("hSparseV2SARec_V2", "hSparseV2SARec_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      histos.add("hpt", "hpt", kTH1F, {configThnAxisPt});
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("CentPercentileMCRecHist", "MC Centrality", kTH1F, {{100, 0.0f, 100.0f}});
      histos.add("h2PhiGen2", "Phi meson gen", kTH2F, {configThnAxisPt, configThnAxisCentrality});
      histos.add("h2PhiRec2", "Phi meson Rec", kTH2F, {configThnAxisPt, configThnAxisCentrality});
      histos.add("hSparseKstarMCGenSA", "hSparseKstarMCGenSA", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
      histos.add("hSparseKstarMCGenCosThetaStar_effy", "hSparseKstarMCGenCosThetaStar_effy", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
      histos.add("hSparseKstarMCRecSA", "hSparseKstarMCRecSA", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
      histos.add("hSparseKstarMCRecCosThetaStar_effy", "hSparseKstarMCRecCosThetaStar_effy", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
    }

    if (doClosureTest) {
      histos.add("closure/KstarSameEventCosThetaStar", "K*0 closure: same-event unlike-sign;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {configThnAxisInvMass, configThnAxisPt, {20, -1.0, 1.0, "cos(#vartheta)"}, configrapAxis, configThnAxisCentrality});
      histos.add("closure/KstarMixedEventCosThetaStar", "K*0 closure: mixed-event unlike-sign;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {configThnAxisInvMass, configThnAxisPt, {20, -1.0, 1.0, "cos(#vartheta)"}, configrapAxis, configThnAxisCentrality});
      histos.add("closure/KstarMCRecTrueCosThetaStar", "K*0 closure: true matched reconstructed;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {configThnAxisInvMass, configThnAxisPt, {20, -1.0, 1.0, "cos(#vartheta)"}, configrapAxis, configThnAxisCentrality});
      histos.add("closure/KstarMCGenCosThetaStar", "K*0 closure: generated;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {configThnAxisInvMass, configThnAxisPt, {20, -1.0, 1.0, "cos(#vartheta)"}, configrapAxis, configThnAxisCentrality});

      const AxisSpec closurePhiMass{phiSA.configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
      const AxisSpec closurePhiPt{phiSA.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec closurePhiCos{20, -1.0, 1.0, "cos(#vartheta)"};
      const AxisSpec closurePhiRap{phiSA.configThnAxisRapidity, "|y|"};
      const AxisSpec closurePhiCent{phiSA.configThnAxisCentrality, "Centrality (%)"};
      histos.add("closure/PhiSameEventCosThetaStar", "phi closure: same-event unlike-sign;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {closurePhiMass, closurePhiPt, closurePhiCos, closurePhiRap, closurePhiCent});
      histos.add("closure/PhiMixedEventCosThetaStar", "phi closure: mixed-event unlike-sign;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {closurePhiMass, closurePhiPt, closurePhiCos, closurePhiRap, closurePhiCent});
      histos.add("closure/PhiMCRecTrueCosThetaStar", "phi closure: true matched reconstructed;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {closurePhiMass, closurePhiPt, closurePhiCos, closurePhiRap, closurePhiCent});
      histos.add("closure/PhiMCGenCosThetaStar", "phi closure: generated;M (GeV/c^{2});p_{T} (GeV/c);cos(#vartheta);|y|;centrality", HistType::kTHnSparseF,
                 {closurePhiMass, closurePhiPt, closurePhiCos, closurePhiRap, closurePhiCent});
    }

    if (doprocessMCkstarWeight) {
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("hImpactParameter", "Impact parameter", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("hEventPlaneAngle", "hEventPlaneAngle", kTH1F, {{200, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}});
      histos.add("hSparseKstarMCGenWeight", "hSparseKstarMCGenWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, configThnAxisPt, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecWeight", "hSparseKstarMCRecWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, configThnAxisPt, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCGenKaonWeight", "hSparseKstarMCGenKaonWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecKaonWeight", "hSparseKstarMCRecKaonWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecKaonMissMatchWeight", "hSparseKstarMCRecKaonMissMatchWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCGenPionWeight", "hSparseKstarMCGenPionWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecPionWeight", "hSparseKstarMCRecPionWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecPionMissMatchWeight", "hSparseKstarMCRecPionMissMatchWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    }

    if (doprocessSE && additionalQAplots1) {
      histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
      histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
      histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
      histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
      histos.add("ResTrackSPFT0CTPC", "ResTrackSPFT0CTPC", kTH2F, {centAxis, resAxis});
      histos.add("ResTrackSPFT0CFT0A", "ResTrackSPFT0CFT0A", kTH2F, {centAxis, resAxis});
      histos.add("ResTrackSPFT0ATPC", "ResTrackSPFT0ATPC", kTH2F, {centAxis, resAxis});
      histos.add("hQFT0CvsCent", "q_{FT0C} vs centrality", kTH2F, {centAxis, qvecAxis});
      histos.add("hQFT0CvsVz", "q_{FT0C} vs Z_{vtx}", kTH2F, {vzAxis, qvecAxis});
      histos.add("hQFT0AvsCent", "q_{FT0A} vs centrality", kTH2F, {centAxis, qvecAxis});
      histos.add("hQFT0AvsVz", "q_{FT0A} vs Z_{vtx}", kTH2F, {vzAxis, qvecAxis});
      histos.add("hQTPCvsCent", "q_{TPC} vs centrality", kTH2F, {centAxis, qvecAxis});
      histos.add("hQTPCvsVz", "q_{TPC} vs Z_{vtx}", kTH2F, {vzAxis, qvecAxis});
    }
    if (doprocessSE && additionalQAplots) {
      // DCA QA
      histos.add("QAbefore/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});

      // PID QA before cuts
      histos.add("QAbefore/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAbefore/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAafter/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAafter/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});

      // DCA QA
      histos.add("QAbefore/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      // PID QA before cuts
      histos.add("QAbefore/TOF_TPC_Mapka_allpi", "TOF + TPC Combined PID for pion;#sigma_{TOF}^{pion};#sigma_{TPC}^{pion}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allpi", "TOF NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAbefore/TPC_Nsigma_allpi", "TPC NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Mapka_allpi", "TOF + TPC Combined PID for pion;#sigma_{TOF}^{pion};#sigma_{TPC}^{pion}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAafter/TOF_Nsigma_allpi", "TOF NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAafter/TPC_Nsigma_allpi", "TPC NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
    }

    if (doprocessSE) {
      histos.add("hMassSameEventLikeNN", "Same-event like-sign (--) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
      histos.add("hMassSameEventLikePP", "Same-event like-sign (++) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
    }
    if (doprocessMixedEvent) {
      histos.add("hMassMixedEventLikeNN", "Mixed-event like-sign (--) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
      histos.add("hMassMixedEventLikePP", "Mixed-event like-sign (++) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
      histos.add("hMassMixedEventUnlike", "Mixed-event unlike-sign mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
    }

    // phi(1020) spin alignment, same names and binning as phipbpb.cxx
    if (doprocessSEPhi || doprocessMEPhi || doprocessMCPhi) {
      const AxisSpec thnAxisInvMassPhi{phiSA.configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPtPhi{phiSA.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec thnAxisCosThetaStarPhi{phiSA.configThnAxisCosThetaStar, "cos(#vartheta_{OP})"};
      const AxisSpec thnAxisCentralityPhi{phiSA.configThnAxisCentrality, "Centrality (%)"};
      const AxisSpec thnAxisRapidityPhi{phiSA.configThnAxisRapidity, "Rapidity"};
      const AxisSpec thnAxisSAPhi{phiSA.configThnAxisSA, "SA"};
      if (cfgSAFrame.value == kRandomEventPlane) {
        histos.add("phi/hPsiRandom", "Random event plane angle", kTH2F, {centAxis, phiAxis});
      }
      if (doprocessSEPhi) {
        histos.add("phi/hSparseV2SameEventSA", "hSparseV2SameEventSA", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisSAPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        histos.add("phi/hSparseV2SameEventCosThetaStar", "hSparseV2SameEventCosThetaStar", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisCosThetaStarPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
      }
      if (doprocessMEPhi) {
        histos.add("phi/hSparseV2MixedEventSA", "hSparseV2MixedEventSA", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisSAPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        histos.add("phi/hSparseV2MixedEventCosThetaStar", "hSparseV2MixedEventCosThetaStar", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisCosThetaStarPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
      }
      if (doprocessMCPhi) {
        histos.add("phi/hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
        histos.add("phi/h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
        histos.add("phi/CentPercentileMCRecHist", "MC Centrality", kTH1F, {{100, 0.0f, 100.0f}});
        histos.add("phi/hSparseV2MCGenSA", "hSparseV2SameEventSA", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisSAPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        histos.add("phi/hSparseV2MCGenCosThetaStar_effy", "hSparseV2SameEventCosThetaStar_effy", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisCosThetaStarPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        histos.add("phi/hSparseV2MCRecSA", "hSparseV2SameEventSA", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisSAPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        histos.add("phi/hSparseV2MCRecCosThetaStar_effy", "hSparseV2SameEventCosThetaStar_effy", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisCosThetaStarPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
      }
    }

    if (doprocessEvtLossSigLossMC) {
      const AxisSpec impactParAxis{evtSigLoss.axisImpactPar, "Impact parameter (fm)"};
      const AxisSpec multMCAxis{evtSigLoss.axisMultMC, "N_{ch}^{gen} (|#eta| < 0.5)"};
      const AxisSpec centLossAxis{100, 0.0, 100.0, "Centrality FT0C (%)"};
      const AxisSpec ptKstarAxis{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec ptPhiAxis{phiSA.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      histos.add("evtSigLoss/MCEventHist", "MC event statistics", kTH1F, {{4, 0.5, 4.5}});
      auto hstat = histos.get<TH1>(HIST("evtSigLoss/MCEventHist"));
      hstat->GetXaxis()->SetBinLabel(1, "All MC events");
      hstat->GetXaxis()->SetBinLabel(2, "MC events with selected reco event");
      hstat->GetXaxis()->SetBinLabel(3, "MC events with no reco event");
      hstat->GetXaxis()->SetBinLabel(4, "MC events with >1 reco event");
      // event loss
      histos.add("evtSigLoss/hImpactParameterGen", "Impact parameter, all generated events", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParameterRec", "Impact parameter, generated events with selected reco event", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParameterGenNoReco", "Impact parameter, generated events with no reco event", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParVsCentRec", "Impact parameter vs centrality of the selected reco event", kTH2F, {centLossAxis, impactParAxis});
      histos.add("evtSigLoss/hMultEta05Gen", "N_{ch}^{gen}, all generated events", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05Rec", "N_{ch}^{gen}, generated events with selected reco event", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05GenNoReco", "N_{ch}^{gen}, generated events with no reco event", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05VsCentRec", "N_{ch}^{gen} vs centrality of the selected reco event", kTH2F, {centLossAxis, multMCAxis});
      // signal loss: K*0 + anti-K*0
      histos.add("evtSigLoss/hKstarGenBeforeEvtSel", "Generated K*0, all events", kTH2F, {ptKstarAxis, impactParAxis});
      histos.add("evtSigLoss/hKstarGenAfterEvtSel", "Generated K*0, events with selected reco event", kTH2F, {ptKstarAxis, impactParAxis});
      histos.add("evtSigLoss/hKstarGenVsMultBeforeEvtSel", "Generated K*0, all events", kTH2F, {ptKstarAxis, multMCAxis});
      histos.add("evtSigLoss/hKstarGenVsMultAfterEvtSel", "Generated K*0, events with selected reco event", kTH2F, {ptKstarAxis, multMCAxis});
      histos.add("evtSigLoss/hKstarGenAfterEvtSelVsCent", "Generated K*0, events with selected reco event", kTH2F, {ptKstarAxis, centLossAxis});
      // signal loss: phi(1020)
      histos.add("evtSigLoss/hPhiGenBeforeEvtSel", "Generated phi, all events", kTH2F, {ptPhiAxis, impactParAxis});
      histos.add("evtSigLoss/hPhiGenAfterEvtSel", "Generated phi, events with selected reco event", kTH2F, {ptPhiAxis, impactParAxis});
      histos.add("evtSigLoss/hPhiGenVsMultBeforeEvtSel", "Generated phi, all events", kTH2F, {ptPhiAxis, multMCAxis});
      histos.add("evtSigLoss/hPhiGenVsMultAfterEvtSel", "Generated phi, events with selected reco event", kTH2F, {ptPhiAxis, multMCAxis});
      histos.add("evtSigLoss/hPhiGenAfterEvtSelVsCent", "Generated phi, events with selected reco event", kTH2F, {ptPhiAxis, centLossAxis});
    }

    // MC closure test of the cos(theta*) chain
    if (doprocessMCClosure || doprocessMCClosureME || doprocessMCClosurePhi || doprocessMCClosureMEPhi) {
      if (closure.rho00Inj.value < 0.f || closure.rho00Inj.value > 1.f) {
        LOGF(fatal, "closure.rho00Inj = %f is outside [0,1]", closure.rho00Inj.value);
      }
      const AxisSpec cosAxisClosure{closure.axisCosThetaStar, "cos(#vartheta*)"};
      const AxisSpec cosFineAxis{100, -1.0, 1.0, "cos(#vartheta*)"};
      if (doprocessMCClosure || doprocessMCClosureME) {
        const AxisSpec massAxis{configThnAxisInvMass, "#it{M}_{K#pi} (GeV/#it{c}^{2})"};
        const AxisSpec ptAxis{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
        const AxisSpec centAxisClosureKstar{configThnAxisCentrality, "Centrality (%)"};
        if (doprocessMCClosure) {
          // bins: 0 all MC events, 1 no reco event, 2 more than one reco event, 3 selected (eff sample), 4 selected (test sample)
          histos.add("closureKstar/hMC", "MC event statistics", kTH1F, {{5, 0.0f, 5.0f}});
          // Acc x Eff sample
          histos.add("closureKstar/hSparseGenEff", "Generated K*, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstar/hSparseRecTrueEff", "Reconstructed true K*, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          // test sample
          histos.add("closureKstar/hSparseGen", "Generated K* (input), test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstar/hSparseRecTrue", "Reconstructed true K*, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstar/hSparseSE", "Same-event unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          if (closure.fillRotBkg) {
            histos.add("closureKstar/hSparseRot", "Rotated unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          }
          histos.add("closureKstar/hCosThetaStarGenVsRec", "cos(#vartheta*) gen vs rec, true K*", kTH2F, {cosFineAxis, cosFineAxis});
        }
        if (doprocessMCClosureME) {
          histos.add("closureKstar/hSparseME", "Mixed-event unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
        }
      }
      if (doprocessMCClosurePhi || doprocessMCClosureMEPhi) {
        const AxisSpec massAxis{phiSA.configThnAxisInvMass, "#it{M}_{KK} (GeV/#it{c}^{2})"};
        const AxisSpec ptAxis{phiSA.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
        const AxisSpec centAxisClosurePhi{phiSA.configThnAxisCentrality, "Centrality (%)"};
        if (doprocessMCClosurePhi) {
          histos.add("closurePhi/hMC", "MC event statistics", kTH1F, {{5, 0.0f, 5.0f}});
          histos.add("closurePhi/hSparseGenEff", "Generated phi, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhi/hSparseRecTrueEff", "Reconstructed true phi, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhi/hSparseGen", "Generated phi (input), test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhi/hSparseRecTrue", "Reconstructed true phi, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhi/hSparseSE", "Same-event K^{+}K^{-} pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          if (closure.fillRotBkg) {
            histos.add("closurePhi/hSparseRot", "Rotated K^{+}K^{-} pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          }
          histos.add("closurePhi/hCosThetaStarGenVsRec", "cos(#vartheta*) gen vs rec, true phi", kTH2F, {cosFineAxis, cosFineAxis});
        }
        if (doprocessMCClosureMEPhi) {
          histos.add("closurePhi/hSparseME", "Mixed-event K^{+}K^{-} pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
        }
      }
    }

    // MC closure test of the cos(theta*) chain, using the reconstructed (finite-resolution) Psi_FT0C instead of the MC true event plane.
    // Needs EPCalibrationTables on the MC sample (i.e. the epvector producer run on it), unlike the true-EP closure above.
    if (doprocessMCClosureRecoEP || doprocessMCClosureMERecoEP || doprocessMCClosureRecoEPPhi || doprocessMCClosureMERecoEPPhi) {
      if (closure.rho00Inj.value < 0.f || closure.rho00Inj.value > 1.f) {
        LOGF(fatal, "closure.rho00Inj = %f is outside [0,1]", closure.rho00Inj.value);
      }
      const AxisSpec cosAxisClosure{closure.axisCosThetaStar, "cos(#vartheta*)"};
      const AxisSpec cosFineAxis{100, -1.0, 1.0, "cos(#vartheta*)"};
      if (doprocessMCClosureRecoEP || doprocessMCClosureMERecoEP) {
        const AxisSpec massAxis{configThnAxisInvMass, "#it{M}_{K#pi} (GeV/#it{c}^{2})"};
        const AxisSpec ptAxis{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
        const AxisSpec centAxisClosureKstar{configThnAxisCentrality, "Centrality (%)"};
        if (doprocessMCClosureRecoEP) {
          // bins: 0 all MC events, 1 no reco event, 2 more than one reco event, 3 selected (eff sample), 4 selected (test sample)
          histos.add("closureKstarRecoEP/hMC", "MC event statistics", kTH1F, {{5, 0.0f, 5.0f}});
          // Acc x Eff sample
          histos.add("closureKstarRecoEP/hSparseGenEff", "Generated K*, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstarRecoEP/hSparseRecTrueEff", "Reconstructed true K*, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          // test sample
          histos.add("closureKstarRecoEP/hSparseGen", "Generated K* (input), test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstarRecoEP/hSparseRecTrue", "Reconstructed true K*, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstarRecoEP/hSparseSE", "Same-event unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          if (closure.fillRotBkg) {
            histos.add("closureKstarRecoEP/hSparseRot", "Rotated unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          }
          histos.add("closureKstarRecoEP/hCosThetaStarGenVsRec", "cos(#vartheta*) gen vs rec, true K*", kTH2F, {cosFineAxis, cosFineAxis});
          // sub-event resolution of the closure sample, needed to derive R2 and un-smear rho00 downstream; same definition as the data ResFT0C*SP
          histos.add("closureKstarRecoEP/ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
          // plain (unweighted) sub-event resolution, as needed by the EP-method correction formula; same definition as the data ResFT0CTPC etc.
          histos.add("closureKstarRecoEP/ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
          // closure-only: Psi_FT0C vs the MC true event plane, to cross-check R from sub-events against the true R
          histos.add("closureKstarRecoEP/ResTrueFT0C", "ResTrueFT0C", kTH2F, {centAxis, resAxis});
        }
        if (doprocessMCClosureMERecoEP) {
          histos.add("closureKstarRecoEP/hSparseME", "Mixed-event unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
        }
      }
      if (doprocessMCClosureRecoEPPhi || doprocessMCClosureMERecoEPPhi) {
        const AxisSpec massAxis{phiSA.configThnAxisInvMass, "#it{M}_{KK} (GeV/#it{c}^{2})"};
        const AxisSpec ptAxis{phiSA.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
        const AxisSpec centAxisClosurePhi{phiSA.configThnAxisCentrality, "Centrality (%)"};
        if (doprocessMCClosureRecoEPPhi) {
          histos.add("closurePhiRecoEP/hMC", "MC event statistics", kTH1F, {{5, 0.0f, 5.0f}});
          histos.add("closurePhiRecoEP/hSparseGenEff", "Generated phi, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhiRecoEP/hSparseRecTrueEff", "Reconstructed true phi, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhiRecoEP/hSparseGen", "Generated phi (input), test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhiRecoEP/hSparseRecTrue", "Reconstructed true phi, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          histos.add("closurePhiRecoEP/hSparseSE", "Same-event K^{+}K^{-} pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          if (closure.fillRotBkg) {
            histos.add("closurePhiRecoEP/hSparseRot", "Rotated K^{+}K^{-} pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
          }
          histos.add("closurePhiRecoEP/hCosThetaStarGenVsRec", "cos(#vartheta*) gen vs rec, true phi", kTH2F, {cosFineAxis, cosFineAxis});
          histos.add("closurePhiRecoEP/ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
          histos.add("closurePhiRecoEP/ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
          histos.add("closurePhiRecoEP/ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
          histos.add("closurePhiRecoEP/ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
          histos.add("closurePhiRecoEP/ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
          histos.add("closurePhiRecoEP/ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
          histos.add("closurePhiRecoEP/ResTrueFT0C", "ResTrueFT0C", kTH2F, {centAxis, resAxis});
        }
        if (doprocessMCClosureMERecoEPPhi) {
          histos.add("closurePhiRecoEP/hSparseME", "Mixed-event K^{+}K^{-} pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosurePhi}, true);
        }
      }
    }

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  static constexpr double MassKa = o2::constants::physics::MassKPlus;
  static constexpr double MassPi = o2::constants::physics::MassPiMinus;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }
  static constexpr float TPCOnlyPt = 0.5f;
  static constexpr int Strategy = 2;
  template <typename T>
  bool selectionPID2(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
        return true;
      }
    }
    if (PID == 1) {
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
        return true;
      }
    }
    return false;
  }
  template <typename T>
  bool strategySelectionPID(const T& candidate, int PID, int strategy)
  {
    if (PID == 0) {
      if (strategy == 0) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == 1) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == Strategy) {
        if (candidate.pt() < TPCOnlyPt && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    if (PID == 1) {
      if (strategy == 0) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == 1) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && ((candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) + (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi())) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == Strategy) {
        if (candidate.pt() < TPCOnlyPt && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    return false;
  }

  double getPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result += o2::constants::math::PI;
    }
    while (result >= o2::constants::math::PI) { // >= not >
      result -= o2::constants::math::PI;
    }
    return result;
  }

  // event-level plane angle for SA: reconstructed Psi2 (event plane) or a random angle (random event plane)
  double getSAEventAngle(double psiEP)
  {
    if (cfgSAFrame.value == kRandomEventPlane) {
      return rndGen.Uniform(-o2::constants::math::PIHalf, o2::constants::math::PIHalf);
    }
    return psiEP;
  }

  // in-plane angle of the SA plane, used for the cos(2(phi* - Psi)) variable
  double getSAPlaneAngle(const ROOT::Math::PxPyPzMVector& mother, double psiSA)
  {
    if (cfgSAFrame.value == kProductionPlane) {
      return mother.Phi(); // production plane is spanned by the beam axis and the mother momentum
    }
    return psiSA;
  }

  // phi(1020) -> K+K- with one kaon taken as the pion: true if the KK mass is inside the veto window (always false if cfgPhiVeto is off)
  bool isPhiReflection(const ROOT::Math::PxPyPzMVector& kaon, const ROOT::Math::PxPyPzMVector& pion)
  {
    if (!cfgPhiVeto) {
      return false;
    }
    ROOT::Math::PxPyPzMVector kaonAsKaon(kaon.Px(), kaon.Py(), kaon.Pz(), MassKa);
    ROOT::Math::PxPyPzMVector pionAsKaon(pion.Px(), pion.Py(), pion.Pz(), MassKa);
    return std::abs((kaonAsKaon + pionAsKaon).M() - o2::constants::physics::MassPhi) < cfgPhiVetoWindow;
  }

  // quantization axis = normal to the SA plane
  ROOT::Math::XYZVector getSAAxis(const ROOT::Math::PxPyPzMVector& mother, double psiSA)
  {
    if (cfgSAFrame.value == kProductionPlane) {
      return ROOT::Math::XYZVector(0., 0., 1.).Cross(mother.Vect()).Unit(); // z x p
    }
    return {std::sin(cfgEPNormalHarmonic.value * psiSA), -std::cos(cfgEPNormalHarmonic.value * psiSA), 0.};
  }

  // ---------------- event selection, common to all K* and phi(1020) process functions ----------------
  template <typename TCollision>
  bool selectionEventRCT(const TCollision& collision)
  {
    return !rctCut.requireRCTFlagChecker || rctCut.rctChecker(collision);
  }

  // additionalEvSel1-4 (data and reconstructed MC)
  template <typename TCollision>
  bool selectionEventBits(const TCollision& collision)
  {
    if (additionalEvSel1 && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (additionalEvSel2 && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (additionalEvSel3 && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    return true;
  }

  // sel8, EP trigger and additionalEvSel1-4
  template <typename TCollision>
  bool selectionEventCuts(const TCollision& collision)
  {
    return collision.sel8() && collision.triggereventep() && selectionEventBits(collision);
  }

  // full data event selection
  template <typename TCollision>
  bool selectionEvent(const TCollision& collision)
  {
    return selectionEventRCT(collision) && selectionEventCuts(collision);
  }

  // mixed-event pair: both events selected and not from the same bunch crossing
  template <typename TCollision>
  bool selectionEventPairME(const TCollision& collision1, const TCollision& collision2)
  {
    if (!selectionEvent(collision1) || !selectionEvent(collision2)) {
      return false;
    }
    return collision1.bcId() != collision2.bcId();
  }

  // ---------------- phi(1020) track and PID selections, identical to phianalysisrun3pbpb.cxx ----------------

  template <typename T>
  bool selectionTrackPhi(const T& candidate)
  {
    if (phiSA.iscustomDCAcut && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > phiSA.cfgITScluster)) {
      return false;
    }
    if (phiSA.ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < phiSA.cfgCutDCAxy && std::abs(candidate.dcaZ()) < phiSA.cfgCutDCAz && candidate.itsNCls() > phiSA.cfgITScluster)) {
      return false;
    }
    return true;
  }

  // kaon PID: circular TPC-TOF cut if TOF is present, TPC-only otherwise (or always, if isNoTOF)
  template <typename T>
  bool selectionPIDPhi(const T& candidate)
  {
    if (!phiSA.isNoTOF && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (phiSA.nsigmaCutCombined * phiSA.nsigmaCutCombined)) {
      return true;
    }
    if (!phiSA.isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < phiSA.nsigmaCutTPC) {
      return true;
    }
    if (phiSA.isNoTOF && std::abs(candidate.tpcNSigmaKa()) < phiSA.nsigmaCutTPC) {
      return true;
    }
    return false;
  }

  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPairPhi(const T1& candidate1, const T2& candidate2)
  {
    double pt1 = candidate1.pt(), pt2 = candidate2.pt();
    double pz1 = candidate1.pz(), pz2 = candidate2.pz();
    double p1 = candidate1.p(), p2 = candidate2.p();
    double angle = std::acos(std::clamp((pt1 * pt2 + pz1 * pz2) / (p1 * p2), -1.0, 1.0)); // clamp = TMath::ACos behaviour
    return !phiSA.isDeepAngle || angle >= phiSA.cfgDeepAngle;
  }

  template <typename T>
  bool isFakeKaonPhi(T const& track)
  {
    return std::abs(track.p() - track.tpcInnerParam()) > phiSA.confFakeKaonCut;
  }

  // phipbpb.cxx builds the EP-frame MC efficiency with the in-plane vector at Psi = 0 (x axis); kept for identical results
  ROOT::Math::XYZVector getSAAxisPhiMC(const ROOT::Math::PxPyPzMVector& mother, double psiSA)
  {
    if (cfgSAFrame.value == kEventPlane) {
      return {std::cos(2.0 * psiSA), std::sin(2.0 * psiSA), 0.};
    }
    return getSAAxis(mother, psiSA);
  }

  template <typename T>
  bool isFakeKaon(T const& track, int /*PID*/)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    return std::abs(pglobal - ptpc) > confFakeKaonCut;
  }

  static constexpr float HalfPI = o2::constants::math::PI * 0.5f;
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle",
                               {6, -HalfPI, HalfPI},
                               "event plane angle"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TH2D* hweight = nullptr;
  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    // scratch vectors: local to this call, never carry state across events/tracks
    ROOT::Math::PxPyPzMVector kstarMother, fourVecDauCM, daughter1, daughter2, kaonrot, kstarrot;
    ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVecNorm;
    ROOT::Math::PxPyPzMVector fourVecDauCMrot;
    ROOT::Math::XYZVector threeVecDauCMrot, threeVecDauCMXYrot;
    double v2 = 0.;
    double v2Rot = 0.;

    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if (!selectionEventRCT(collision)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 1.5);
    if (!selectionEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 2.5);
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto qFT0C = collision.qFT0C();
    auto qFT0A = collision.qFT0A();
    auto qTPC = collision.qTPC();
    histos.fill(HIST("hEvtSelInfo"), 3.5);
    if (additionalQAplots1) {
      histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
      histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
      histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
      histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
      histos.fill(HIST("ResFT0CTPC"), centrality, std::cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0A"), centrality, std::cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("hVtxZ"), collision.posZ());
      histos.fill(HIST("hQFT0CvsCent"), centrality, qFT0C);
      histos.fill(HIST("hQFT0CvsVz"), collision.posZ(), qFT0C);
      histos.fill(HIST("hQFT0AvsCent"), centrality, qFT0A);
      histos.fill(HIST("hQFT0AvsVz"), collision.posZ(), qFT0A);
      histos.fill(HIST("hQTPCvsCent"), centrality, qTPC);
      histos.fill(HIST("hQTPCvsVz"), collision.posZ(), qTPC);
    }
    // one angle per event: Psi2 from FT0C, or a random angle for the random-EP null test
    auto psiSA = getSAEventAngle(psiFT0C);
    if (fillSA && cfgSAFrame.value == kRandomEventPlane) {
      histos.fill(HIST("hPsiRandom"), centrality, psiSA);
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    if (useWeight && (currentRunNumber != lastRunNumber)) {
      hweight = ccdb->getForTimeStamp<TH2D>(confWeightPath.value, bc.timestamp());
    }
    lastRunNumber = currentRunNumber;
    float weight1 = 1.0;
    float weight2 = 1.0;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }

      auto track1ID = track1.globalIndex();
      if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
        continue;
      }
      if (isTOFOnly && !selectionPID2(track1, 0)) {
        continue;
      }

      if (useWeight) {
        if (track1.pt() < cfgMaxTrackPt &&
            track1.pt() > cfgMinTrackPt) {
          weight1 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track1.pt() + 0.000005)) * std::cos(2.0 * getPhiInRange(track1.phi() - psiFT0C));
        } else {
          weight1 = 1;
        }
      }
      for (const auto& track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }

        auto track2ID = track2.globalIndex();
        if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track2, 1)) {
          continue;
        }

        if (track2ID == track1ID) {
          continue;
        }

        if (additionalQAplots) {
          histos.fill(HIST("QAafter/TPC_Nsigma_allka"), track1.pt(), track1.tpcNSigmaKa(), centrality);
          histos.fill(HIST("QAafter/TOF_Nsigma_allka"), track1.pt(), track1.tofNSigmaKa(), centrality);
          histos.fill(HIST("QAafter/trkDCAxyka"), track1.dcaXY());
          histos.fill(HIST("QAafter/trkDCAzka"), track1.dcaZ());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_allka"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_allpi"), track2.tofNSigmaPi(), track2.tpcNSigmaPi());
          histos.fill(HIST("QAafter/TPC_Nsigma_allpi"), track2.pt(), track2.tpcNSigmaPi(), centrality);
          histos.fill(HIST("QAafter/TOF_Nsigma_allpi"), track2.pt(), track2.tofNSigmaPi(), centrality);
          histos.fill(HIST("QAafter/trkDCAxypi"), track2.dcaXY());
          histos.fill(HIST("QAafter/trkDCAzpi"), track2.dcaZ());
        }
        if (useWeight) {
          if (track2.pt() < cfgMaxTrackPt &&
              track2.pt() > cfgMinTrackPt) {
            weight2 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track2.pt() + 0.000005)) * std::cos(2.0 * getPhiInRange(track2.phi() - psiFT0C));
          } else {
            weight2 = 1;
          }
        }
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
        kstarMother = daughter1 + daughter2;
        if (std::abs(kstarMother.Rapidity()) > confRapidity) {
          continue;
        }
        if (isPhiReflection(daughter1, daughter2)) { // phi(1020) veto (cfgPhiVeto); also skips the rotations of this pair
          continue;
        }
        auto phiMinusPsi = getPhiInRange(kstarMother.Phi() - psiFT0C);

        if (useSP) {
          v2 = std::cos(2.0 * phiMinusPsi) * qFT0C;
        }
        if (!useSP) {
          v2 = std::cos(2.0 * phiMinusPsi);
        }
        auto totalweight = weight1 * weight2;
        static constexpr float MinTotalWeight = 5e-7f;
        if (totalweight <= MinTotalWeight) {
          totalweight = 1.0;
        }
        if (additionalQAplots1) {
          histos.fill(HIST("ResTrackSPFT0CTPC"), centrality, qFT0C * qTPC * std::cos(2.0 * (psiFT0C - psiTPC)));
          histos.fill(HIST("ResTrackSPFT0CFT0A"), centrality, qFT0C * qFT0A * std::cos(2.0 * (psiFT0C - psiFT0A)));
          histos.fill(HIST("ResTrackSPFT0ATPC"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        }
        if (!fillSA) {

          if (useWeight) {
            histos.fill(HIST("hSparseV2SASameEvent_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality, 1 / totalweight);
          } else {
            histos.fill(HIST("hSparseV2SASameEvent_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
          }
        }
        int track1Sign = track1.sign();
        int track2Sign = track2.sign();

        if (track1Sign * track2Sign > 0) {
          if (track1Sign > 0) {
            histos.fill(HIST("hMassSameEventLikePP"), kstarMother.M(), centrality);
          } else {
            histos.fill(HIST("hMassSameEventLikeNN"), kstarMother.M(), centrality);
          }
        } else {
        }

        if (fillSA) {
          ROOT::Math::Boost boost{kstarMother.BoostToCM()};
          fourVecDauCM = boost(daughter1);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          eventplaneVecNorm = getSAAxis(kstarMother, psiSA);
          auto cosPhistarminuspsi = getPhiInRange(fourVecDauCM.Phi() - getSAPlaneAngle(kstarMother, psiSA));
          auto sa = std::cos(2.0 * cosPhistarminuspsi);
          auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());

          if (doClosureTest && track1Sign * track2Sign < 0) {
            histos.fill(HIST("closure/KstarSameEventCosThetaStar"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, std::abs(kstarMother.Rapidity()), centrality);
          }

          if (track1Sign * track2Sign < 0) {
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapsameunlike"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, kstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapsameunlike"), kstarMother.M(), kstarMother.Pt(), sa, kstarMother.Rapidity(), centrality);
            }
          } else if (fillLikeSignSA && track1Sign * track2Sign > 0) {
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapsamelike"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, kstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapsamelike"), kstarMother.M(), kstarMother.Pt(), sa, kstarMother.Rapidity(), centrality);
            }
          }
        }
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            if (!fillSA) {
              histos.fill(HIST("hRotation"), rotangle);
            }
            auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
            auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
            kaonrot = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), MassKa);
            kstarrot = kaonrot + daughter2;
            if (std::abs(kstarrot.Rapidity()) > confRapidity) {
              continue;
            }
            if (isPhiReflection(kaonrot, daughter2)) { // phi(1020) veto (cfgPhiVeto) on the rotated pair
              continue;
            }
            auto phiMinusPsiRot = getPhiInRange(kstarrot.Phi() - psiFT0C);

            if (useSP) {
              v2Rot = std::cos(2.0 * phiMinusPsiRot) * qFT0C;
            }
            if (!useSP) {
              v2Rot = std::cos(2.0 * phiMinusPsiRot);
            }
            if (!fillSA) {
              histos.fill(HIST("hSparseV2SASameEventRotational_V2"), kstarrot.M(), kstarrot.Pt(), v2Rot, centrality);
            }
            if (fillSA) {
              if (track1Sign * track2Sign < 0) {
                ROOT::Math::Boost boost{kstarrot.BoostToCM()};
                fourVecDauCMrot = boost(kaonrot);
                threeVecDauCMrot = fourVecDauCMrot.Vect();
                threeVecDauCMXYrot = ROOT::Math::XYZVector(threeVecDauCMrot.X(), threeVecDauCMrot.Y(), 0.);
                // production-plane axis must follow the rotated candidate
                auto saAxisRot = getSAAxis(kstarrot, psiSA);
                auto cosPhistarminuspsirot = getPhiInRange(fourVecDauCMrot.Phi() - getSAPlaneAngle(kstarrot, psiSA));
                auto sarot = std::cos(2.0 * cosPhistarminuspsirot);
                auto cosThetaStarrot = saAxisRot.Dot(threeVecDauCMrot) / std::sqrt(threeVecDauCMrot.Mag2()) / std::sqrt(saAxisRot.Mag2());
                if (usepolar) {
                  histos.fill(HIST("hSparseSAvsraprot"), kstarrot.M(), kstarrot.Pt(), cosThetaStarrot, kstarrot.Rapidity(), centrality);
                } else {
                  histos.fill(HIST("hSparseSAvsraprot"), kstarrot.M(), kstarrot.Pt(), sarot, kstarrot.Rapidity(), centrality);
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processSE, "Process Same event latest", true);

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    // scratch vectors: local to this call, never carry state across pairs
    ROOT::Math::PxPyPzMVector kstarMother, fourVecDauCM, daughter1, daughter2;
    ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVecNorm;
    double v2 = 0.;

    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!selectionEventPairME(collision1, collision2)) {
        continue;
      }
      auto centrality = collision1.centFT0C();
      auto psiFT0C1 = collision1.psiFT0C();
      auto qFT0C1 = collision1.qFT0C();
      auto psiFT0C2 = collision2.psiFT0C();
      auto qFT0C2 = collision2.qFT0C();
      auto psiSA1 = getSAEventAngle(psiFT0C1);

      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track1, 0)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track2, 1)) {
          continue;
        }

        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);

        kstarMother = daughter1 + daughter2;
        if (std::abs(kstarMother.Rapidity()) > confRapidity) {
          continue;
        }
        if (isPhiReflection(daughter1, daughter2)) { // phi(1020) veto (cfgPhiVeto), same as in the same event
          continue;
        }

        int s1 = track1.sign();
        int s2 = track2.sign();

        if (s1 * s2 < 0) {

          auto phi1 = track1.phi();
          auto phi2 = track2.phi();
          auto phiKstar = kstarMother.Phi();

          double term1 = qFT0C1 * std::cos(2.0 * getPhiInRange(phi1 - psiFT0C1)) * std::cos(2.0 * getPhiInRange(phi1 - phiKstar));
          double term2 = qFT0C2 * std::cos(2.0 * getPhiInRange(phi2 - psiFT0C2)) * std::cos(2.0 * getPhiInRange(phi2 - phiKstar));

          v2 = term1 + term2;

          if (!fillSA) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
          }
          histos.fill(HIST("hMassMixedEventUnlike"), kstarMother.M(), centrality);

          if (fillSA) {
            ROOT::Math::Boost boost{kstarMother.BoostToCM()};
            fourVecDauCM = boost(daughter1);
            threeVecDauCM = fourVecDauCM.Vect();
            threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
            eventplaneVecNorm = getSAAxis(kstarMother, psiSA1);
            auto cosPhistarminuspsi = getPhiInRange(fourVecDauCM.Phi() - getSAPlaneAngle(kstarMother, psiSA1));
            auto sa = std::cos(2.0 * cosPhistarminuspsi);
            auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
            if (doClosureTest) {
              histos.fill(HIST("closure/KstarMixedEventCosThetaStar"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, std::abs(kstarMother.Rapidity()), centrality);
            }
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapmix"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, kstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapmix"), kstarMother.M(), kstarMother.Pt(), sa, kstarMother.Rapidity(), centrality);
            }
          }
        } else {

          if (s1 > 0) {
            histos.fill(HIST("hMassMixedEventLikePP"), kstarMother.M(), centrality);
          } else {
            histos.fill(HIST("hMassMixedEventLikeNN"), kstarMother.M(), centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMixedEvent, "Process Mixed event", true);

  void processMC(CollisionMCTrueTable::iterator const& /*TrueCollision*/, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    // scratch vectors: local to this call, never carry state across particles
    ROOT::Math::PxPyPzMVector kstarMother, kaonPlus, pionMinus;
    double v2 = 0.;

    histos.fill(HIST("hMC"), 0);
    if (RecCollisions.size() == 0) {
      histos.fill(HIST("hMC"), 1);
      return;
    }
    if (RecCollisions.size() > 1) {
      histos.fill(HIST("hMC"), 2);
      return;
    }
    for (const auto& RecCollision : RecCollisions) {
      auto psiFT0C = 0.0;
      histos.fill(HIST("hMC"), 3);
      if (!RecCollision.sel8()) {
        histos.fill(HIST("hMC"), 4);
        continue;
      }
      if (!selectionEventBits(RecCollision)) {
        continue;
      }
      histos.fill(HIST("hMC"), 5);
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("hMC"), 6);
        continue;
      }
      histos.fill(HIST("hMC"), 7);
      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("CentPercentileMCRecHist"), centrality);
      auto psiSA = getSAEventAngle(psiFT0C); // same angle for rec and gen of this event
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (const auto& track1 : rectrackspart) {
        if (!selectionTrack(track1)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track1, 0)) {
          continue;
        }
        if (!track1.has_mcParticle()) {
          continue;
        }
        auto track1ID = track1.index();
        for (const auto& track2 : rectrackspart) {
          auto track2ID = track2.index();
          // roles are fixed (track1 = kaon candidate, track2 = pion candidate, checked against the MC PDG below),
          // so every true K pi pair is found exactly once whatever the track order; only skip the same track
          if (track2ID == track1ID) {
            continue;
          }
          if (!selectionTrack(track2)) {
            continue;
          }
          if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
            continue;
          }
          if (isTOFOnly && !selectionPID2(track2, 1)) {
            continue;
          }
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (track1.sign() * track2.sign() > 0) {
            continue;
          }
          const auto mctrack1 = track1.mcParticle();
          const auto mctrack2 = track2.mcParticle();
          int track1PDG = std::abs(mctrack1.pdgCode());
          int track2PDG = std::abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (track1PDG != PDG_t::kKPlus || track2PDG != PDG_t::kPiPlus) {
            continue;
          }
          for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (std::abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              // K*0 and anti-K*0, as in data (all unlike-sign K pi pairs)
              if (pdgcheck && std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kK0Star892) {
                continue;
              }
              if (!isTOFOnly && !(strategySelectionPID(track1, 0, strategyPID) || strategySelectionPID(track2, 1, strategyPID))) {
                continue;
              }
              if (isTOFOnly && !(selectionPID2(track1, 0) || selectionPID2(track2, 1))) {
                continue;
              }
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              // track1 is always the kaon-PID candidate and track2 always the pion-PID candidate (enforced above);
              // assign unconditionally so the boost below always uses the true kaon, regardless of its charge sign
              kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
              pionMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
              kstarMother = kaonPlus + pionMinus;
              if (std::abs(kstarMother.Rapidity()) > confRapidity) {
                continue;
              }
              if (isPhiReflection(kaonPlus, pionMinus)) { // phi(1020) veto (cfgPhiVeto), so the efficiency includes the true K* it removes
                continue;
              }
              auto phiMinusPsi = getPhiInRange(kstarMother.Phi() - psiFT0C);

              v2 = std::cos(2.0 * phiMinusPsi);

              histos.fill(HIST("hSparseV2SARec_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
              histos.fill(HIST("h2PhiRec2"), kstarMother.pt(), centrality);
              histos.fill(HIST("hpt"), kstarMother.Pt());

              {
                ROOT::Math::Boost boost{kstarMother.BoostToCM()};
                auto fourVecDauCMRec = boost(kaonPlus);
                auto threeVecDauCMRec = fourVecDauCMRec.Vect();
                auto eventplaneVecNormRec = getSAAxis(kstarMother, psiSA);
                auto cosPhistarminuspsiRec = getPhiInRange(fourVecDauCMRec.Phi() - getSAPlaneAngle(kstarMother, psiSA));
                auto saRec = std::cos(2.0 * cosPhistarminuspsiRec);
                auto cosThetaStarRec = eventplaneVecNormRec.Dot(threeVecDauCMRec) / std::sqrt(threeVecDauCMRec.Mag2()) / std::sqrt(eventplaneVecNormRec.Mag2());

                histos.fill(HIST("hSparseKstarMCRecSA"), kstarMother.M(), kstarMother.Pt(), saRec, std::abs(kstarMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseKstarMCRecCosThetaStar_effy"), kstarMother.M(), kstarMother.Pt(), cosThetaStarRec, std::abs(kstarMother.Rapidity()), centrality);
                if (doClosureTest) {
                  histos.fill(HIST("closure/KstarMCRecTrueCosThetaStar"), kstarMother.M(), kstarMother.Pt(), cosThetaStarRec, std::abs(kstarMother.Rapidity()), centrality);
                }
              }
            }
          }
        }
      }
      // loop over generated particle
      for (const auto& mcParticle : GenParticles) {
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        // K*0 and anti-K*0, as in data and in the reconstructed loop above
        if (pdgcheck && std::abs(mcParticle.pdgCode()) != o2::constants::physics::kK0Star892) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        static constexpr std::size_t NumberOfDaughters = 2;

        if (kDaughters.size() != NumberOfDaughters) {
          continue;
        }
        // daughtp = kaon daughter found, daughtm = pion daughter found (either charge: K+ pi- or K- pi+);
        // kaonPlus / pionMinus hold the kaon / pion daughter, so the boost uses the kaon as in the reconstructed loop
        auto daughtp = false;
        auto daughtm = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kKPlus) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            if (!genacceptancecut) {
              daughtp = true;
            }
            kaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), MassKa);
          } else if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kPiPlus) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            if (!genacceptancecut) {
              daughtm = true;
            }
            pionMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), MassPi);
          }
        }
        if (daughtp && daughtm) {
          kstarMother = kaonPlus + pionMinus;
          if (std::abs(kstarMother.Rapidity()) > confRapidity) {
            continue;
          }
          auto phiMinusPsi = getPhiInRange(kstarMother.Phi() - psiFT0C);

          v2 = std::cos(2.0 * phiMinusPsi);

          histos.fill(HIST("hSparseV2SAGen_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
          histos.fill(HIST("h2PhiGen2"), kstarMother.pt(), centrality);

          {
            ROOT::Math::Boost boost{kstarMother.BoostToCM()};
            auto fourVecDauCMGen = boost(kaonPlus);
            auto threeVecDauCMGen = fourVecDauCMGen.Vect();
            auto eventplaneVecNormGen = getSAAxis(kstarMother, psiSA);
            auto cosPhistarminuspsiGen = getPhiInRange(fourVecDauCMGen.Phi() - getSAPlaneAngle(kstarMother, psiSA));
            auto saGen = std::cos(2.0 * cosPhistarminuspsiGen);
            auto cosThetaStarGen = eventplaneVecNormGen.Dot(threeVecDauCMGen) / std::sqrt(threeVecDauCMGen.Mag2()) / std::sqrt(eventplaneVecNormGen.Mag2());

            histos.fill(HIST("hSparseKstarMCGenSA"), kstarMother.M(), kstarMother.Pt(), saGen, std::abs(kstarMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseKstarMCGenCosThetaStar_effy"), kstarMother.M(), kstarMother.Pt(), cosThetaStarGen, std::abs(kstarMother.Rapidity()), centrality);
            if (doClosureTest) {
              histos.fill(HIST("closure/KstarMCGenCosThetaStar"), kstarMother.M(), kstarMother.Pt(), cosThetaStarGen, std::abs(kstarMother.Rapidity()), centrality);
            }
          }
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(Kstarpbpb, processMC, "Process MC", false);

  void processMCkstarWeight(CollisionMCTrueTable::iterator const& TrueCollision, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    // scratch vectors: local to this call, never carry state across particles
    ROOT::Math::PxPyPzMVector kaonPlus, pionMinus;

    float imp = TrueCollision.impactParameter();
    float evPhi = TrueCollision.eventPlaneAngle() / 2.0;
    static constexpr std::array<float, 10> CentEdges = {
      0.0f, 3.49f, 4.93f, 6.98f, 8.55f,
      9.87f, 11.0f, 12.1f, 13.1f, 14.0f};

    static constexpr std::array<float, 9> CentValues = {
      2.5f, 7.5f, 15.0f, 25.0f, 35.0f,
      45.0f, 55.0f, 65.0f, 75.0f};
    float centclass = -999.f;

    for (size_t i = 0; i < CentValues.size(); ++i) {
      if (imp >= CentEdges[i] && imp < CentEdges[i + 1]) {
        centclass = CentValues[i];
        break;
      }
    }
    histos.fill(HIST("hImpactParameter"), imp);
    histos.fill(HIST("hEventPlaneAngle"), evPhi);
    static constexpr float MinCentrality = 0.0f;
    static constexpr float MaxCentrality = 80.0f;

    if (centclass < MinCentrality || centclass > MaxCentrality) {
      return;
    }
    for (const auto& RecCollision : RecCollisions) {
      auto psiFT0C = TrueCollision.eventPlaneAngle();
      /*
  if (!RecCollision.sel8()) {
        continue;
      }
      if (!RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      continue;
      }
      */
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        continue;
      }
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (const auto& track1 : rectrackspart) {
        if (!track1.has_mcParticle()) {
          continue;
        }

        const auto mctrack1 = track1.mcParticle();

        if (selectionTrack(track1) && strategySelectionPID(track1, 0, strategyPID) && std::abs(mctrack1.pdgCode()) == PDG_t::kKPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecKaonWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && track1.pt() > TPCOnlyPt && track1.hasTOF() && std::abs(track1.tofNSigmaKa()) > nsigmaCutTOF && std::abs(track1.tpcNSigmaKa()) < nsigmaCutTPC && std::abs(mctrack1.pdgCode()) == PDG_t::kKPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecKaonMissMatchWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && strategySelectionPID(track1, 1, strategyPID) && std::abs(mctrack1.pdgCode()) == PDG_t::kPiPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecPionWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && track1.pt() > TPCOnlyPt && track1.hasTOF() && std::abs(track1.tofNSigmaPi()) > nsigmaCutTOF && std::abs(track1.tpcNSigmaPi()) < nsigmaCutTPC && std::abs(mctrack1.pdgCode()) == PDG_t::kPiPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecPionMissMatchWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        auto track1ID = track1.index();
        for (const auto& track2 : rectrackspart) {
          if (!track2.has_mcParticle()) {
            continue;
          }
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          const auto mctrack2 = track2.mcParticle();
          int track1PDG = std::abs(mctrack1.pdgCode());
          int track2PDG = std::abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (track1PDG != PDG_t::kKPlus || track2PDG != PDG_t::kPiPlus) {
            continue;
          }
          for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (std::abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kK0Star892) {
                continue;
              }
              // if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              if (avoidsplitrackMC && oldindex == mothertrack1.index()) {
                histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              // oldindex = mothertrack1.globalIndex();
              oldindex = mothertrack1.index();
              auto phiMinusPsi = getPhiInRange(mothertrack1.phi() - psiFT0C);
              histos.fill(HIST("hSparseKstarMCRecWeight"), centclass, phiMinusPsi, std::pow(std::cos(2.0 * phiMinusPsi), 2.0), mothertrack1.pt(), mothertrack1.eta());
            }
          }
        }
      }
      // loop over generated particle
      for (const auto& mcParticle : GenParticles) {
        static constexpr float MaxEtaAcceptance = 0.8f;
        if (std::abs(mcParticle.eta()) > MaxEtaAcceptance) {
          continue;
        }
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus && mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCGenKaonWeight"), centclass, getPhiInRange(mcParticle.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mcParticle.phi() - psiFT0C)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCGenPionWeight"), centclass, getPhiInRange(mcParticle.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mcParticle.phi() - psiFT0C)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (mcParticle.pdgCode() != o2::constants::physics::kK0Star892) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        static constexpr std::size_t NumberOfDaughters = 2;

        if (kDaughters.size() != NumberOfDaughters) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == +PDG_t::kKPlus) {
            if (kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            kaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), MassKa);
          } else if (kCurrentDaughter.pdgCode() == -PDG_t::kPiPlus) {
            if (kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            pionMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), MassPi);
          }
        }
        if (daughtp && daughtm) {
          auto phiMinusPsiGen = getPhiInRange(mcParticle.phi() - psiFT0C);
          histos.fill(HIST("hSparseKstarMCGenWeight"), centclass, phiMinusPsiGen, std::pow(std::cos(2.0 * phiMinusPsiGen), 2.0), mcParticle.pt(), mcParticle.eta());
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(Kstarpbpb, processMCkstarWeight, "Process MC kstar Weight", false);

  // ================= phi(1020) spin alignment: port of the SA part of phipbpb.cxx =================
  struct SAValues {
    double cosThetaStar;
    double sa;
  };

  // K- in the phi rest frame (as in phipbpb.cxx) w.r.t. the chosen quantization axis
  SAValues getSAValuesPhi(const ROOT::Math::PxPyPzMVector& mother, const ROOT::Math::PxPyPzMVector& kaonMinus, const ROOT::Math::XYZVector& axis, double psiSA)
  {
    ROOT::Math::Boost boost{mother.BoostToCM()};
    auto threeVecDau = boost(kaonMinus).Vect();
    auto cosThetaStar = axis.Dot(threeVecDau) / std::sqrt(threeVecDau.Mag2()) / std::sqrt(axis.Mag2());
    auto sa = std::cos(2.0 * getPhiInRange(threeVecDau.Phi() - getSAPlaneAngle(mother, psiSA)));
    return {.cosThetaStar = cosThetaStar, .sa = sa};
  }

  void processSEPhi(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    if (!selectionEvent(collision)) {
      return;
    }
    auto centrality = collision.centFT0C();
    auto psiSA = getSAEventAngle(collision.psiFT0C());
    if (cfgSAFrame.value == kRandomEventPlane) {
      histos.fill(HIST("phi/hPsiRandom"), centrality, psiSA);
    }

    ROOT::Math::PxPyPzMVector kaonPlusPhi, kaonMinusPhi, phiMother;
    for (const auto& track1 : tracks) {
      if (!(track1.signed1Pt() > cfgCutCharge.value)) { // positive kaon
        continue;
      }
      if (!selectionTrackPhi(track1) || !selectionPIDPhi(track1)) {
        continue;
      }
      for (const auto& track2 : tracks) {
        if (!(track2.signed1Pt() < cfgCutCharge.value)) { // negative kaon
          continue;
        }
        if (!selectionTrackPhi(track2) || !selectionPIDPhi(track2)) {
          continue;
        }
        if (!selectionPairPhi(track1, track2)) {
          continue;
        }
        if (phiSA.removeFakeTrack && (isFakeKaonPhi(track1) || isFakeKaonPhi(track2))) {
          continue;
        }
        kaonPlusPhi = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
        kaonMinusPhi = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
        phiMother = kaonPlusPhi + kaonMinusPhi;
        auto absRapidity = std::abs(phiMother.Rapidity());
        if (absRapidity > phiSA.confRapidity) {
          continue;
        }
        auto [cosThetaStar, sa] = getSAValuesPhi(phiMother, kaonMinusPhi, getSAAxis(phiMother, psiSA), psiSA);
        if (doClosureTest) {
          histos.fill(HIST("closure/PhiSameEventCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, absRapidity, centrality);
        }
        histos.fill(HIST("phi/hSparseV2SameEventSA"), phiMother.M(), phiMother.Pt(), sa, absRapidity, centrality);
        histos.fill(HIST("phi/hSparseV2SameEventCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, absRapidity, centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processSEPhi, "Process same event phi(1020) spin alignment", false);

  void processMEPhi(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    ROOT::Math::PxPyPzMVector kaonPlusPhi, kaonMinusPhi, phiMother;
    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!selectionEventPairME(collision1, collision2)) {
        continue;
      }
      auto centrality = collision1.centFT0C();
      auto psiSA = getSAEventAngle(collision1.psiFT0C());
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrackPhi(track1) || !selectionTrackPhi(track2)) {
          continue;
        }
        if (!selectionPIDPhi(track1) || !selectionPIDPhi(track2)) {
          continue;
        }
        if (!selectionPairPhi(track1, track2)) {
          continue;
        }
        if (phiSA.removeFakeTrack && (isFakeKaonPhi(track1) || isFakeKaonPhi(track2))) {
          continue;
        }
        if (track1.sign() > 0) {
          kaonPlusPhi = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonMinusPhi = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
        } else {
          kaonMinusPhi = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonPlusPhi = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
        }
        phiMother = kaonPlusPhi + kaonMinusPhi;
        auto absRapidity = std::abs(phiMother.Rapidity());
        if (absRapidity > phiSA.confRapidity) {
          continue;
        }
        auto [cosThetaStar, sa] = getSAValuesPhi(phiMother, kaonMinusPhi, getSAAxis(phiMother, psiSA), psiSA);
        if (doClosureTest) {
          histos.fill(HIST("closure/PhiMixedEventCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, absRapidity, centrality);
        }
        histos.fill(HIST("phi/hSparseV2MixedEventSA"), phiMother.M(), phiMother.Pt(), sa, absRapidity, centrality);
        histos.fill(HIST("phi/hSparseV2MixedEventCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, absRapidity, centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMEPhi, "Process mixed event phi(1020) spin alignment", false);

  void processMCPhi(CollisionMCTrueTable::iterator const& /*TrueCollision*/, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    histos.fill(HIST("phi/hMC"), 0);
    if (RecCollisions.size() == 0) {
      histos.fill(HIST("phi/hMC"), 1);
      return;
    }
    if (RecCollisions.size() > 1) {
      histos.fill(HIST("phi/hMC"), 2);
      return;
    }
    ROOT::Math::PxPyPzMVector kaonPlusPhi, kaonMinusPhi, phiMother;
    for (const auto& RecCollision : RecCollisions) {
      auto psiFT0C = 0.0;
      histos.fill(HIST("phi/hMC"), 3);
      if (!RecCollision.sel8()) {
        histos.fill(HIST("phi/hMC"), 4);
        continue;
      }
      if (!selectionEventBits(RecCollision)) {
        continue;
      }
      histos.fill(HIST("phi/hMC"), 5); // same bin meaning as the K* hMC
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("phi/hMC"), 6);
        continue;
      }
      histos.fill(HIST("phi/hMC"), 7);
      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("phi/CentPercentileMCRecHist"), centrality);
      auto psiSA = getSAEventAngle(psiFT0C); // same angle for rec and gen of this event
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (const auto& track1 : rectrackspart) {
        if (!selectionTrackPhi(track1)) {
          continue;
        }
        if (!selectionPIDPhi(track1)) {
          continue;
        }
        if (!track1.has_mcParticle()) {
          continue;
        }
        auto track1ID = track1.index();
        for (const auto& track2 : rectrackspart) {
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          if (!selectionTrackPhi(track2)) {
            continue;
          }
          if (!selectionPIDPhi(track2)) {
            continue;
          }
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (!selectionPairPhi(track1, track2)) {
            continue;
          }
          if (track1.sign() * track2.sign() > 0) {
            continue;
          }
          const auto mctrack1 = track1.mcParticle();
          const auto mctrack2 = track2.mcParticle();
          if (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (std::abs(mctrack1.pdgCode()) != PDG_t::kKPlus || std::abs(mctrack2.pdgCode()) != PDG_t::kKPlus) {
            continue;
          }
          for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (std::abs(mothertrack1.y()) > phiSA.confRapidity) {
                continue;
              }
              if (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kPhi) {
                continue;
              }
              if (phiSA.avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                histos.fill(HIST("phi/h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              if (track1.sign() > 0) {
                kaonPlusPhi = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
                kaonMinusPhi = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
              } else {
                kaonMinusPhi = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
                kaonPlusPhi = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
              }
              phiMother = kaonPlusPhi + kaonMinusPhi;
              // reconstructed-pair rapidity, same cut as in processSEPhi / processMEPhi (true-mother y is cut above)
              if (std::abs(phiMother.Rapidity()) > phiSA.confRapidity) {
                continue;
              }
              auto [cosThetaStar, sa] = getSAValuesPhi(phiMother, kaonMinusPhi, getSAAxisPhiMC(phiMother, psiSA), psiSA);
              histos.fill(HIST("phi/hSparseV2MCRecCosThetaStar_effy"), phiMother.M(), phiMother.Pt(), cosThetaStar, std::abs(phiMother.Rapidity()), centrality);
              if (doClosureTest) {
                histos.fill(HIST("closure/PhiMCRecTrueCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, std::abs(phiMother.Rapidity()), centrality);
              }
              histos.fill(HIST("phi/hSparseV2MCRecSA"), phiMother.M(), phiMother.Pt(), sa, std::abs(phiMother.Rapidity()), centrality);
            }
          }
        }
      }
      // loop over generated particle
      for (const auto& mcParticle : GenParticles) {
        if (std::abs(mcParticle.y()) > phiSA.confRapidity) {
          continue;
        }
        if (mcParticle.pdgCode() != o2::constants::physics::kPhi) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        static constexpr std::size_t NumberOfDaughters = 2;
        if (kDaughters.size() != NumberOfDaughters) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          bool inAcceptance = !phiSA.genacceptancecut || (kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta);
          if (kCurrentDaughter.pdgCode() == +PDG_t::kKPlus) {
            daughtp = daughtp || inAcceptance;
            kaonPlusPhi = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), MassKa);
          } else if (kCurrentDaughter.pdgCode() == -PDG_t::kKPlus) {
            daughtm = daughtm || inAcceptance;
            kaonMinusPhi = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), MassKa);
          }
        }
        if (daughtp && daughtm) {
          phiMother = kaonPlusPhi + kaonMinusPhi;
          auto [cosThetaStar, sa] = getSAValuesPhi(phiMother, kaonMinusPhi, getSAAxisPhiMC(phiMother, psiSA), psiSA);
          histos.fill(HIST("phi/hSparseV2MCGenCosThetaStar_effy"), phiMother.M(), phiMother.Pt(), cosThetaStar, std::abs(phiMother.Rapidity()), centrality);
          if (doClosureTest) {
            histos.fill(HIST("closure/PhiMCGenCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, std::abs(phiMother.Rapidity()), centrality);
          }
          histos.fill(HIST("phi/hSparseV2MCGenSA"), phiMother.M(), phiMother.Pt(), sa, std::abs(phiMother.Rapidity()), centrality);
        }
      }
    } // rec collision loop
  }
  PROCESS_SWITCH(Kstarpbpb, processMCPhi, "Process MC phi(1020) spin alignment", false);

  // ================= event loss and signal loss for K* and phi(1020) =================
  // adapted from processEvtLossSigLossMC1 in phianalysisrun3pbpb.cxx. The accepted MC event uses the same reconstructed-event
  // selection as processMC / processMCPhi (sel8, additionalEvSel1-4, |vz|) and, by default, the same single-reco-collision
  // requirement; the generated K* and phi use the same rapidity and daughter conditions as the generated loops there.
  void processEvtLossSigLossMC(McCollisionMults::iterator const& mcCollision, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles)
  {
    if (evtSigLoss.cutVzGen && std::abs(mcCollision.posZ()) > cfgCutVertex) {
      return;
    }
    if (evtSigLoss.isApplyInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    if (evtSigLoss.isApplyTVX && (mcCollision.multMCFT0C() <= 0 || mcCollision.multMCFT0A() <= 0)) {
      return;
    }
    const float impactPar = mcCollision.impactParameter();
    const float multMC = mcCollision.multMCNParticlesEta05();

    // all generated events
    histos.fill(HIST("evtSigLoss/MCEventHist"), 1);
    histos.fill(HIST("evtSigLoss/hImpactParameterGen"), impactPar);
    histos.fill(HIST("evtSigLoss/hMultEta05Gen"), multMC);
    if (RecCollisions.size() == 0) {
      histos.fill(HIST("evtSigLoss/MCEventHist"), 3);
      histos.fill(HIST("evtSigLoss/hImpactParameterGenNoReco"), impactPar);
      histos.fill(HIST("evtSigLoss/hMultEta05GenNoReco"), multMC);
    }
    if (RecCollisions.size() > 1) {
      histos.fill(HIST("evtSigLoss/MCEventHist"), 4);
    }

    // generated events with a selected reconstructed collision (event loss)
    bool selected = false;
    float centrality = -999.f;
    if (!evtSigLoss.requireSingleReco || RecCollisions.size() == 1) {
      for (const auto& RecCollision : RecCollisions) {
        if (!RecCollision.sel8() || !selectionEventBits(RecCollision) || std::abs(RecCollision.posZ()) > cfgCutVertex) {
          continue;
        }
        selected = true;
        centrality = RecCollision.centFT0C();
      }
    }
    if (selected) {
      histos.fill(HIST("evtSigLoss/MCEventHist"), 2);
      histos.fill(HIST("evtSigLoss/hImpactParameterRec"), impactPar);
      histos.fill(HIST("evtSigLoss/hMultEta05Rec"), multMC);
      histos.fill(HIST("evtSigLoss/hImpactParVsCentRec"), centrality, impactPar);
      histos.fill(HIST("evtSigLoss/hMultEta05VsCentRec"), centrality, multMC);
    }

    // generated K*0 / anti-K*0 and phi(1020) (signal loss)
    static constexpr std::size_t NumberOfDaughters = 2;
    for (const auto& mcParticle : GenParticles) {
      const int pdgMother = std::abs(mcParticle.pdgCode());
      if (pdgMother != o2::constants::physics::kK0Star892 && mcParticle.pdgCode() != o2::constants::physics::kPhi) {
        continue;
      }
      const bool isKstar = (pdgMother == o2::constants::physics::kK0Star892);
      if (std::abs(mcParticle.y()) > (isKstar ? confRapidity.value : phiSA.confRapidity.value)) {
        continue;
      }
      auto daughters = mcParticle.daughters_as<aod::McParticles>();
      if (daughters.size() != NumberOfDaughters) {
        continue;
      }
      bool hasKaon = false, hasPion = false, hasKPlus = false, hasKMinus = false;
      ROOT::Math::PxPyPzMVector dau1, dau2;
      for (const auto& dau : daughters) {
        if (!dau.isPhysicalPrimary()) {
          continue;
        }
        if (isKstar) {
          if (std::abs(dau.pdgCode()) == PDG_t::kKPlus) {
            hasKaon = true;
            dau1 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          } else if (std::abs(dau.pdgCode()) == PDG_t::kPiPlus) {
            hasPion = true;
            dau2 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassPi);
          }
        } else {
          if (dau.pdgCode() == PDG_t::kKPlus) {
            hasKPlus = true;
            dau1 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          } else if (dau.pdgCode() == PDG_t::kKMinus) {
            hasKMinus = true;
            dau2 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          }
        }
      }
      const double pt = (dau1 + dau2).Pt();
      if (isKstar && hasKaon && hasPion) {
        histos.fill(HIST("evtSigLoss/hKstarGenBeforeEvtSel"), pt, impactPar);
        histos.fill(HIST("evtSigLoss/hKstarGenVsMultBeforeEvtSel"), pt, multMC);
        if (selected) {
          histos.fill(HIST("evtSigLoss/hKstarGenAfterEvtSel"), pt, impactPar);
          histos.fill(HIST("evtSigLoss/hKstarGenVsMultAfterEvtSel"), pt, multMC);
          histos.fill(HIST("evtSigLoss/hKstarGenAfterEvtSelVsCent"), pt, centrality);
        }
      } else if (!isKstar && hasKPlus && hasKMinus) {
        histos.fill(HIST("evtSigLoss/hPhiGenBeforeEvtSel"), pt, impactPar);
        histos.fill(HIST("evtSigLoss/hPhiGenVsMultBeforeEvtSel"), pt, multMC);
        if (selected) {
          histos.fill(HIST("evtSigLoss/hPhiGenAfterEvtSel"), pt, impactPar);
          histos.fill(HIST("evtSigLoss/hPhiGenVsMultAfterEvtSel"), pt, multMC);
          histos.fill(HIST("evtSigLoss/hPhiGenAfterEvtSelVsCent"), pt, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processEvtLossSigLossMC, "Process event loss and signal loss for K* and phi(1020)", false);

  // ================= MC closure test of the cos(theta*) chain =================
  // Test sample (reconstructed MC analysed like data):
  //   hSparseSE      all same-event unlike-sign pairs (signal + background), true pairs weighted by W
  //   hSparseME      mixed-event unlike-sign pairs (weight 1)                -> step 1: SE - norm * ME
  //   hSparseRecTrue reconstructed true K* (phi), weighted by W              -> step 2: truth for the raw yield
  //   hSparseGen     generated K* (phi), daughters in acceptance, weighted   -> step 4: input yield and rho00
  // Efficiency sample (independent events if closure.splitSample):
  //   hSparseRecTrueEff / hSparseGenEff = Acc x Eff(cos theta*)              -> step 3
  // W(cos) = 3/2 [(1 - rho00) + (3 rho00 - 1) cos^2] with the GENERATED cos(theta*) (W = 1 if injectPolarization is off)

  // W(cos) = 3/2 [(1 - rho) + (3 rho - 1) cos^2]; <W> = 1 over cos in [-1,1]; W = 1 for rho = 1/3
  double getInjectionWeight(double cosThetaStar)
  {
    if (!closure.injectPolarization) {
      return 1.0;
    }
    const double rho = closure.rho00Inj.value;
    return 1.5 * ((1.0 - rho) + (3.0 * rho - 1.0) * cosThetaStar * cosThetaStar);
  }

  // cos(theta*) of 'dau' in the rest frame of 'mother' with respect to 'axis'
  double getCosThetaStar(const ROOT::Math::PxPyPzMVector& mother, const ROOT::Math::PxPyPzMVector& dau, const ROOT::Math::XYZVector& axis)
  {
    ROOT::Math::Boost boost{mother.BoostToCM()};
    auto threeVecDau = boost(dau).Vect();
    return axis.Dot(threeVecDau) / std::sqrt(threeVecDau.Mag2()) / std::sqrt(axis.Mag2());
  }

  // one Psi per MC event (call once per event: the random-EP option draws a new angle at each call)
  double getClosurePsi(double trueEventPlaneAngle)
  {
    return getSAEventAngle(closure.useTrueEP ? trueEventPlaneAngle : 0.0);
  }

  // even MC events -> Acc x Eff sample, odd MC events -> test sample (both if splitSample is off)
  static constexpr int64_t SplitSampleModulus = 2;
  bool isClosureEffSample(int64_t mcCollisionIndex)
  {
    return !closure.splitSample || (mcCollisionIndex % SplitSampleModulus == 0);
  }
  bool isClosureTestSample(int64_t mcCollisionIndex)
  {
    return !closure.splitSample || (mcCollisionIndex % SplitSampleModulus == 1);
  }

  // track enters the SE / ME / rotated pairs: always if requirePrimaryTracks is off, else only if matched to a physical primary
  template <typename TTrack>
  bool selectionClosurePrimary(const TTrack& track)
  {
    return !closure.requirePrimaryTracks || (track.has_mcParticle() && track.mcParticle().isPhysicalPrimary());
  }

  // number of reco collisions per MC collision (for the single-reco requirement of the ME, same as the SE)
  std::vector<int> countRecoPerMcCollision(aod::McCollisionLabels const& mcLabels, int64_t nMcCollisions)
  {
    std::vector<int> nRecPerMc(nMcCollisions, 0);
    for (const auto& label : mcLabels) {
      if (label.mcCollisionId() >= 0 && label.mcCollisionId() < nMcCollisions) {
        nRecPerMc[label.mcCollisionId()]++;
      }
    }
    return nRecPerMc;
  }

  // mixed-event pair of the closure test sample: same event class as the SE
  template <typename TCollision>
  bool selectionClosurePairME(const TCollision& collision1, const TCollision& collision2, const std::vector<int>& nRecPerMc)
  {
    if (!selectionEventMC(collision1) || !selectionEventMC(collision2) || collision1.bcId() == collision2.bcId()) {
      return false;
    }
    if (!collision1.has_mcCollision() || !collision2.has_mcCollision()) {
      return false;
    }
    if (!isClosureTestSample(collision1.mcCollisionId()) || !isClosureTestSample(collision2.mcCollisionId())) {
      return false;
    }
    if (closure.meRequireSingleReco && (nRecPerMc[collision1.mcCollisionId()] != 1 || nRecPerMc[collision2.mcCollisionId()] != 1)) {
      return false;
    }
    return true;
  }

  // reco-event selection of processMC / processMCPhi, plus the data centrality cut (identical for SE and ME)
  template <typename TCollision>
  bool selectionEventMC(const TCollision& collision)
  {
    return collision.sel8() && selectionEventBits(collision) && std::abs(collision.posZ()) < cfgCutVertex && std::abs(collision.centFT0C()) < cfgCutCentrality;
  }

  // reco-event selection of the reconstructed-Psi_FT0C closure: selectionEventMC plus the EP-calibration trigger required
  // on data (processSE / processSEPhi via selectionEventCuts), only available once EPCalibrationTables is joined
  template <typename TCollision>
  bool selectionEventMCRecoEP(const TCollision& collision)
  {
    return selectionEventMC(collision) && collision.triggereventep();
  }

  // mixed-event pair of the reconstructed-Psi_FT0C closure test sample: same as selectionClosurePairME, with the EP trigger
  template <typename TCollision>
  bool selectionClosurePairMERecoEP(const TCollision& collision1, const TCollision& collision2, const std::vector<int>& nRecPerMc)
  {
    if (!selectionEventMCRecoEP(collision1) || !selectionEventMCRecoEP(collision2) || collision1.bcId() == collision2.bcId()) {
      return false;
    }
    if (!collision1.has_mcCollision() || !collision2.has_mcCollision()) {
      return false;
    }
    if (!isClosureTestSample(collision1.mcCollisionId()) || !isClosureTestSample(collision2.mcCollisionId())) {
      return false;
    }
    if (closure.meRequireSingleReco && (nRecPerMc[collision1.mcCollisionId()] != 1 || nRecPerMc[collision2.mcCollisionId()] != 1)) {
      return false;
    }
    return true;
  }

  // K pi pair from the decay of one primary K*0 / anti-K*0 (|y| < confRapidity)
  template <typename TTrack>
  bool isTrueKstarPair(const TTrack& trackKaon, const TTrack& trackPion)
  {
    if (!trackKaon.has_mcParticle() || !trackPion.has_mcParticle()) {
      return false;
    }
    const auto mcKaon = trackKaon.mcParticle();
    const auto mcPion = trackPion.mcParticle();
    if (!mcKaon.isPhysicalPrimary() || !mcPion.isPhysicalPrimary()) {
      return false;
    }
    if (std::abs(mcKaon.pdgCode()) != PDG_t::kKPlus || std::abs(mcPion.pdgCode()) != PDG_t::kPiPlus) {
      return false;
    }
    for (const auto& motherKaon : mcKaon.template mothers_as<aod::McParticles>()) {
      for (const auto& motherPion : mcPion.template mothers_as<aod::McParticles>()) {
        if (motherKaon != motherPion) {
          continue;
        }
        if (std::abs(motherKaon.pdgCode()) != o2::constants::physics::kK0Star892) {
          continue;
        }
        if (std::abs(motherKaon.y()) > confRapidity) {
          continue;
        }
        return true;
      }
    }
    return false;
  }

  // K+ K- pair from the decay of one primary phi(1020) (|y| < phiSA.confRapidity)
  template <typename TTrack>
  bool isTruePhiPair(const TTrack& trackPlus, const TTrack& trackMinus)
  {
    if (!trackPlus.has_mcParticle() || !trackMinus.has_mcParticle()) {
      return false;
    }
    const auto mcPlus = trackPlus.mcParticle();
    const auto mcMinus = trackMinus.mcParticle();
    if (!mcPlus.isPhysicalPrimary() || !mcMinus.isPhysicalPrimary()) {
      return false;
    }
    if (mcPlus.pdgCode() != PDG_t::kKPlus || mcMinus.pdgCode() != PDG_t::kKMinus) {
      return false;
    }
    for (const auto& motherPlus : mcPlus.template mothers_as<aod::McParticles>()) {
      for (const auto& motherMinus : mcMinus.template mothers_as<aod::McParticles>()) {
        if (motherPlus != motherMinus) {
          continue;
        }
        if (motherPlus.pdgCode() != o2::constants::physics::kPhi) {
          continue;
        }
        if (std::abs(motherPlus.y()) > phiSA.confRapidity) {
          continue;
        }
        return true;
      }
    }
    return false;
  }

  // ---------------- K*0: same event, true reco, generated ----------------
  // track1 = kaon candidate, track2 = pion candidate, unlike-sign, same cuts as processSE.
  // The axis is built with getSAAxis, as in processMC (efficiency) and processSE (data).
  void processMCClosure(CollisionMCTrueTable::iterator const& trueCollision, CollisionMCRecTableCentFT0C const& recCollisions, TrackMCTrueTable const& genParticles, FilTrackMCRecTable const& recTracks)
  {
    ROOT::Math::PxPyPzMVector kaon, pion, kstar, kaonRot, kstarRot, kaonGen, pionGen, kstarGen;
    static constexpr std::size_t NumberOfDaughters = 2;

    histos.fill(HIST("closureKstar/hMC"), 0.5);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("closureKstar/hMC"), 1.5);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("closureKstar/hMC"), 2.5);
      return;
    }
    const bool isEffSample = isClosureEffSample(trueCollision.globalIndex());
    const bool isTestSample = isClosureTestSample(trueCollision.globalIndex());
    const double psi = getClosurePsi(trueCollision.eventPlaneAngle()); // same angle for the rec and the gen part of this event
    for (const auto& recCollision : recCollisions) {
      if (!selectionEventMC(recCollision)) {
        continue;
      }
      if (isEffSample) {
        histos.fill(HIST("closureKstar/hMC"), 3.5);
      }
      if (isTestSample) {
        histos.fill(HIST("closureKstar/hMC"), 4.5);
      }
      const auto centrality = recCollision.centFT0C();
      auto recTracksThisColl = recTracks.sliceBy(perCollision, recCollision.globalIndex());

      // ---- reconstructed pairs ----
      for (const auto& track1 : recTracksThisColl) {
        if (!selectionTrack(track1) || !selectionClosurePrimary(track1)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track1, 0)) {
          continue;
        }
        for (const auto& track2 : recTracksThisColl) {
          if (track2.globalIndex() == track1.globalIndex()) {
            continue;
          }
          if (track1.sign() * track2.sign() > 0) {
            continue;
          }
          if (!selectionTrack(track2) || !selectionClosurePrimary(track2)) {
            continue;
          }
          if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
            continue;
          }
          if (isTOFOnly && !selectionPID2(track2, 1)) {
            continue;
          }
          kaon = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          pion = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
          kstar = kaon + pion;
          if (std::abs(kstar.Rapidity()) > confRapidity) {
            continue;
          }
          if (isPhiReflection(kaon, pion)) {
            continue;
          }
          const double cosRec = getCosThetaStar(kstar, kaon, getSAAxis(kstar, psi));

          double weight = 1.0;
          if (isTrueKstarPair(track1, track2)) {
            const auto mcKaon = track1.mcParticle();
            const auto mcPion = track2.mcParticle();
            kaonGen = ROOT::Math::PxPyPzMVector(mcKaon.px(), mcKaon.py(), mcKaon.pz(), MassKa);
            pionGen = ROOT::Math::PxPyPzMVector(mcPion.px(), mcPion.py(), mcPion.pz(), MassPi);
            kstarGen = kaonGen + pionGen;
            const double cosGen = getCosThetaStar(kstarGen, kaonGen, getSAAxis(kstarGen, psi));
            weight = getInjectionWeight(cosGen);
            // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
            if (isEffSample) {
              histos.fill(HIST("closureKstar/hSparseRecTrueEff"), kstar.M(), kstar.Pt(), cosRec, centrality, 1.0);
            }
            if (isTestSample) {
              histos.fill(HIST("closureKstar/hSparseRecTrue"), kstar.M(), kstar.Pt(), cosRec, centrality, weight);
              histos.fill(HIST("closureKstar/hCosThetaStarGenVsRec"), cosGen, cosRec);
            }
          }
          if (!isTestSample) {
            continue;
          }
          histos.fill(HIST("closureKstar/hSparseSE"), kstar.M(), kstar.Pt(), cosRec, centrality, weight);

          // rotated background, same rotation as in processSE (weight 1)
          if (closure.fillRotBkg && nBkgRotations.value > 1) {
            const double angleStart = confMinRot.value;
            const double angleStep = (confMaxRot.value - angleStart) / (nBkgRotations.value - 1.0);
            for (int irot = 0; irot < nBkgRotations.value; irot++) {
              const double rotAngle = angleStart + irot * angleStep;
              kaonRot = ROOT::Math::PxPyPzMVector(track1.px() * std::cos(rotAngle) - track1.py() * std::sin(rotAngle), track1.px() * std::sin(rotAngle) + track1.py() * std::cos(rotAngle), track1.pz(), MassKa);
              kstarRot = kaonRot + pion;
              if (std::abs(kstarRot.Rapidity()) > confRapidity) {
                continue;
              }
              if (isPhiReflection(kaonRot, pion)) {
                continue;
              }
              const double cosRot = getCosThetaStar(kstarRot, kaonRot, getSAAxis(kstarRot, psi));
              histos.fill(HIST("closureKstar/hSparseRot"), kstarRot.M(), kstarRot.Pt(), cosRot, centrality);
            }
          }
        }
      }

      // ---- generated K*0 / anti-K*0: Acc x Eff denominator (eff sample) and input truth (test sample) ----
      for (const auto& mcParticle : genParticles) {
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (std::abs(mcParticle.pdgCode()) != o2::constants::physics::kK0Star892) {
          continue;
        }
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        if (daughters.size() != NumberOfDaughters) {
          continue;
        }
        bool hasKaon = false;
        bool hasPion = false;
        for (const auto& dau : daughters) {
          if (!dau.isPhysicalPrimary()) {
            continue;
          }
          const bool inAcceptance = !closure.genAcceptanceCut || (dau.pt() > cfgCutPT && std::abs(dau.eta()) < cfgCutEta);
          if (std::abs(dau.pdgCode()) == PDG_t::kKPlus) {
            hasKaon = hasKaon || inAcceptance;
            kaonGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          } else if (std::abs(dau.pdgCode()) == PDG_t::kPiPlus) {
            hasPion = hasPion || inAcceptance;
            pionGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassPi);
          }
        }
        if (!hasKaon || !hasPion) {
          continue;
        }
        kstarGen = kaonGen + pionGen;
        if (std::abs(kstarGen.Rapidity()) > confRapidity) {
          continue;
        }
        const double cosGen = getCosThetaStar(kstarGen, kaonGen, getSAAxis(kstarGen, psi));
        // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
        if (isEffSample) {
          histos.fill(HIST("closureKstar/hSparseGenEff"), kstarGen.M(), kstarGen.Pt(), cosGen, centrality, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closureKstar/hSparseGen"), kstarGen.M(), kstarGen.Pt(), cosGen, centrality, weight);
        }
      }
    } // rec collision loop
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosure, "Process MC closure of cos(theta*) for K*: SE, true rec, gen", false);

  // ---------------- K*0: mixed event of the test sample ----------------
  // Same pairing as processMixedEvent (kaon from event 1, pion from event 2, unlike-sign), axis from event 1.
  void processMCClosureME(CollisionMCRecMixTable const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    ROOT::Math::PxPyPzMVector kaon, pion, kstar;
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeMCClosure binningOnPositions{{axisVertex, axisMultiplicityClass}, true};
    SameKindPair<CollisionMCRecMixTable, FilTrackMCRecTable, BinningTypeMCClosure> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!selectionClosurePairME(collision1, collision2, nRecPerMc)) {
        continue;
      }
      const auto centrality = collision1.centFT0C();
      const double psi = getClosurePsi(collision1.mcCollision().eventPlaneAngle());
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        if (!selectionClosurePrimary(track1) || !selectionClosurePrimary(track2)) {
          continue;
        }
        if (!isTOFOnly && (!strategySelectionPID(track1, 0, strategyPID) || !strategySelectionPID(track2, 1, strategyPID))) {
          continue;
        }
        if (isTOFOnly && (!selectionPID2(track1, 0) || !selectionPID2(track2, 1))) {
          continue;
        }
        kaon = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
        pion = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
        kstar = kaon + pion;
        if (std::abs(kstar.Rapidity()) > confRapidity) {
          continue;
        }
        if (isPhiReflection(kaon, pion)) {
          continue;
        }
        const double cosRec = getCosThetaStar(kstar, kaon, getSAAxis(kstar, psi));
        histos.fill(HIST("closureKstar/hSparseME"), kstar.M(), kstar.Pt(), cosRec, centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureME, "Process MC closure of cos(theta*) for K*: mixed event", false);

  // ---------------- phi(1020): same event, true reco, generated ----------------
  // track1 = K+, track2 = K-, as in processSEPhi; cos(theta*) of the K-.
  // The axis is built with getSAAxis, as in processSEPhi/processMEPhi (the real data analysis).
  void processMCClosurePhi(CollisionMCTrueTable::iterator const& trueCollision, CollisionMCRecTableCentFT0C const& recCollisions, TrackMCTrueTable const& genParticles, FilTrackMCRecTable const& recTracks)
  {
    ROOT::Math::PxPyPzMVector kaonPlus, kaonMinus, phi, kaonPlusRot, phiRot, kaonPlusGen, kaonMinusGen, phiGen;
    static constexpr std::size_t NumberOfDaughters = 2;

    histos.fill(HIST("closurePhi/hMC"), 0.5);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("closurePhi/hMC"), 1.5);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("closurePhi/hMC"), 2.5);
      return;
    }
    const bool isEffSample = isClosureEffSample(trueCollision.globalIndex());
    const bool isTestSample = isClosureTestSample(trueCollision.globalIndex());
    const double psi = getClosurePsi(trueCollision.eventPlaneAngle());
    for (const auto& recCollision : recCollisions) {
      if (!selectionEventMC(recCollision)) {
        continue;
      }
      if (isEffSample) {
        histos.fill(HIST("closurePhi/hMC"), 3.5);
      }
      if (isTestSample) {
        histos.fill(HIST("closurePhi/hMC"), 4.5);
      }
      const auto centrality = recCollision.centFT0C();
      auto recTracksThisColl = recTracks.sliceBy(perCollision, recCollision.globalIndex());

      // ---- reconstructed pairs ----
      for (const auto& track1 : recTracksThisColl) {
        if (!(track1.signed1Pt() > cfgCutCharge.value)) { // K+
          continue;
        }
        if (!selectionTrackPhi(track1) || !selectionPIDPhi(track1) || !selectionClosurePrimary(track1)) {
          continue;
        }
        for (const auto& track2 : recTracksThisColl) {
          if (!(track2.signed1Pt() < cfgCutCharge.value)) { // K-
            continue;
          }
          if (!selectionTrackPhi(track2) || !selectionPIDPhi(track2) || !selectionClosurePrimary(track2)) {
            continue;
          }
          if (!selectionPairPhi(track1, track2)) {
            continue;
          }
          if (phiSA.removeFakeTrack && (isFakeKaonPhi(track1) || isFakeKaonPhi(track2))) {
            continue;
          }
          kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
          phi = kaonPlus + kaonMinus;
          if (std::abs(phi.Rapidity()) > phiSA.confRapidity) {
            continue;
          }
          const double cosRec = getCosThetaStar(phi, kaonMinus, getSAAxis(phi, psi));

          double weight = 1.0;
          if (isTruePhiPair(track1, track2)) {
            const auto mcPlus = track1.mcParticle();
            const auto mcMinus = track2.mcParticle();
            kaonPlusGen = ROOT::Math::PxPyPzMVector(mcPlus.px(), mcPlus.py(), mcPlus.pz(), MassKa);
            kaonMinusGen = ROOT::Math::PxPyPzMVector(mcMinus.px(), mcMinus.py(), mcMinus.pz(), MassKa);
            phiGen = kaonPlusGen + kaonMinusGen;
            const double cosGen = getCosThetaStar(phiGen, kaonMinusGen, getSAAxis(phiGen, psi));
            weight = getInjectionWeight(cosGen);
            // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
            if (isEffSample) {
              histos.fill(HIST("closurePhi/hSparseRecTrueEff"), phi.M(), phi.Pt(), cosRec, centrality, 1.0);
            }
            if (isTestSample) {
              histos.fill(HIST("closurePhi/hSparseRecTrue"), phi.M(), phi.Pt(), cosRec, centrality, weight);
              histos.fill(HIST("closurePhi/hCosThetaStarGenVsRec"), cosGen, cosRec);
            }
          }
          if (!isTestSample) {
            continue;
          }
          histos.fill(HIST("closurePhi/hSparseSE"), phi.M(), phi.Pt(), cosRec, centrality, weight);

          // rotated background: rotate the K+ (weight 1)
          if (closure.fillRotBkg && nBkgRotations.value > 1) {
            const double angleStart = confMinRot.value;
            const double angleStep = (confMaxRot.value - angleStart) / (nBkgRotations.value - 1.0);
            for (int irot = 0; irot < nBkgRotations.value; irot++) {
              const double rotAngle = angleStart + irot * angleStep;
              kaonPlusRot = ROOT::Math::PxPyPzMVector(track1.px() * std::cos(rotAngle) - track1.py() * std::sin(rotAngle), track1.px() * std::sin(rotAngle) + track1.py() * std::cos(rotAngle), track1.pz(), MassKa);
              phiRot = kaonPlusRot + kaonMinus;
              if (std::abs(phiRot.Rapidity()) > phiSA.confRapidity) {
                continue;
              }
              const double cosRot = getCosThetaStar(phiRot, kaonMinus, getSAAxis(phiRot, psi));
              histos.fill(HIST("closurePhi/hSparseRot"), phiRot.M(), phiRot.Pt(), cosRot, centrality);
            }
          }
        }
      }

      // ---- generated phi: Acc x Eff denominator (eff sample) and input truth (test sample) ----
      for (const auto& mcParticle : genParticles) {
        if (std::abs(mcParticle.y()) > phiSA.confRapidity) {
          continue;
        }
        if (mcParticle.pdgCode() != o2::constants::physics::kPhi) {
          continue;
        }
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        if (daughters.size() != NumberOfDaughters) {
          continue;
        }
        bool hasKPlus = false;
        bool hasKMinus = false;
        for (const auto& dau : daughters) {
          if (!dau.isPhysicalPrimary()) {
            continue;
          }
          const bool inAcceptance = !closure.genAcceptanceCut || (dau.pt() > cfgCutPT && std::abs(dau.eta()) < cfgCutEta);
          if (dau.pdgCode() == PDG_t::kKPlus) {
            hasKPlus = hasKPlus || inAcceptance;
            kaonPlusGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          } else if (dau.pdgCode() == PDG_t::kKMinus) {
            hasKMinus = hasKMinus || inAcceptance;
            kaonMinusGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          }
        }
        if (!hasKPlus || !hasKMinus) {
          continue;
        }
        phiGen = kaonPlusGen + kaonMinusGen;
        const double cosGen = getCosThetaStar(phiGen, kaonMinusGen, getSAAxis(phiGen, psi));
        // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
        if (isEffSample) {
          histos.fill(HIST("closurePhi/hSparseGenEff"), phiGen.M(), phiGen.Pt(), cosGen, centrality, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closurePhi/hSparseGen"), phiGen.M(), phiGen.Pt(), cosGen, centrality, weight);
        }
      }
    } // rec collision loop
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosurePhi, "Process MC closure of cos(theta*) for phi(1020): SE, true rec, gen", false);

  // ---------------- phi(1020): mixed event of the test sample ----------------
  // Same pairing as processMEPhi (K+K- from two events), cos(theta*) of the K-, axis from event 1.
  void processMCClosureMEPhi(CollisionMCRecMixTable const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    ROOT::Math::PxPyPzMVector kaonPlus, kaonMinus, phi;
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeMCClosure binningOnPositions{{axisVertex, axisMultiplicityClass}, true};
    SameKindPair<CollisionMCRecMixTable, FilTrackMCRecTable, BinningTypeMCClosure> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!selectionClosurePairME(collision1, collision2, nRecPerMc)) {
        continue;
      }
      const auto centrality = collision1.centFT0C();
      const double psi = getClosurePsi(collision1.mcCollision().eventPlaneAngle());
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrackPhi(track1) || !selectionTrackPhi(track2)) {
          continue;
        }
        if (!selectionClosurePrimary(track1) || !selectionClosurePrimary(track2)) {
          continue;
        }
        if (!selectionPIDPhi(track1) || !selectionPIDPhi(track2)) {
          continue;
        }
        if (!selectionPairPhi(track1, track2)) {
          continue;
        }
        if (phiSA.removeFakeTrack && (isFakeKaonPhi(track1) || isFakeKaonPhi(track2))) {
          continue;
        }
        if (track1.sign() > 0) {
          kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
        } else {
          kaonMinus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonPlus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
        }
        phi = kaonPlus + kaonMinus;
        if (std::abs(phi.Rapidity()) > phiSA.confRapidity) {
          continue;
        }
        const double cosRec = getCosThetaStar(phi, kaonMinus, getSAAxis(phi, psi));
        histos.fill(HIST("closurePhi/hSparseME"), phi.M(), phi.Pt(), cosRec, centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureMEPhi, "Process MC closure of cos(theta*) for phi(1020): mixed event", false);

  // ================= MC closure test of the cos(theta*) chain, reconstructed (finite-resolution) Psi_FT0C =================
  // Same as processMCClosure(Phi)/processMCClosureME(Phi) above, but Psi is the reconstructed Psi_FT0C (as on data) instead of
  // the MC true event plane, and the sub-event resolution of the sample is filled too (ResFT0C*SP, needed to derive R2 and
  // un-smear the fitted raw rho00 downstream). Kept as separate process functions (own histogram directories, own collision
  // tables joined with EPCalibrationTables) so the true-EP closure above keeps zero dependency on the EP-calibration producer.

  // ---------------- K*0: same event, true reco, generated (reconstructed Psi_FT0C) ----------------
  void processMCClosureRecoEP(CollisionMCTrueTable::iterator const& trueCollision, CollisionMCRecTableClosure const& recCollisions, TrackMCTrueTable const& genParticles, FilTrackMCRecTable const& recTracks)
  {
    ROOT::Math::PxPyPzMVector kaon, pion, kstar, kaonRot, kstarRot, kaonGen, pionGen, kstarGen;
    static constexpr std::size_t NumberOfDaughters = 2;

    histos.fill(HIST("closureKstarRecoEP/hMC"), 0.5);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("closureKstarRecoEP/hMC"), 1.5);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("closureKstarRecoEP/hMC"), 2.5);
      return;
    }
    const bool isEffSample = isClosureEffSample(trueCollision.globalIndex());
    const bool isTestSample = isClosureTestSample(trueCollision.globalIndex());
    // true axis: defines the injected weight and the Gen/GenEff truth, so the test is not tautological against the smeared axis
    const double psiTrue = getSAEventAngle(trueCollision.eventPlaneAngle());
    for (const auto& recCollision : recCollisions) {
      if (!selectionEventMCRecoEP(recCollision)) {
        continue;
      }
      if (isEffSample) {
        histos.fill(HIST("closureKstarRecoEP/hMC"), 3.5);
      }
      if (isTestSample) {
        histos.fill(HIST("closureKstarRecoEP/hMC"), 4.5);
      }
      const auto centrality = recCollision.centFT0C();
      const double recoPsiFT0C = recCollision.psiFT0C();
      const double psi = getSAEventAngle(recoPsiFT0C); // measured axis: defines cosRec (SE/ME/RecTrue/RecTrueEff/Rot) only
      {
        const double psiTPC = recCollision.psiTPC();
        const double psiFT0A = recCollision.psiFT0A();
        const double qFT0C = recCollision.qFT0C();
        const double qTPC = recCollision.qTPC();
        const double qFT0A = recCollision.qFT0A();
        histos.fill(HIST("closureKstarRecoEP/ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        // plain (unweighted) sub-event resolution, as needed by the EP-method correction formula; same definition as the data ResFT0CTPC etc.
        histos.fill(HIST("closureKstarRecoEP/ResFT0CTPC"), centrality, std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0CFT0A"), centrality, std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
        // closure-only: Psi_FT0C vs the MC true event plane, to cross-check R from sub-events against the true R
        histos.fill(HIST("closureKstarRecoEP/ResTrueFT0C"), centrality, std::cos(2.0 * (recoPsiFT0C - trueCollision.eventPlaneAngle())));
      }
      auto recTracksThisColl = recTracks.sliceBy(perCollision, recCollision.globalIndex());

      // ---- reconstructed pairs ----
      for (const auto& track1 : recTracksThisColl) {
        if (!selectionTrack(track1) || !selectionClosurePrimary(track1)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track1, 0)) {
          continue;
        }
        for (const auto& track2 : recTracksThisColl) {
          if (track2.globalIndex() == track1.globalIndex()) {
            continue;
          }
          if (track1.sign() * track2.sign() > 0) {
            continue;
          }
          if (!selectionTrack(track2) || !selectionClosurePrimary(track2)) {
            continue;
          }
          if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
            continue;
          }
          if (isTOFOnly && !selectionPID2(track2, 1)) {
            continue;
          }
          kaon = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          pion = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
          kstar = kaon + pion;
          if (std::abs(kstar.Rapidity()) > confRapidity) {
            continue;
          }
          if (isPhiReflection(kaon, pion)) {
            continue;
          }
          const double cosRec = getCosThetaStar(kstar, kaon, getSAAxis(kstar, psi));

          double weight = 1.0;
          if (isTrueKstarPair(track1, track2)) {
            const auto mcKaon = track1.mcParticle();
            const auto mcPion = track2.mcParticle();
            kaonGen = ROOT::Math::PxPyPzMVector(mcKaon.px(), mcKaon.py(), mcKaon.pz(), MassKa);
            pionGen = ROOT::Math::PxPyPzMVector(mcPion.px(), mcPion.py(), mcPion.pz(), MassPi);
            kstarGen = kaonGen + pionGen;
            const double cosGen = getCosThetaStar(kstarGen, kaonGen, getSAAxis(kstarGen, psiTrue));
            weight = getInjectionWeight(cosGen);
            // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
            if (isEffSample) {
              histos.fill(HIST("closureKstarRecoEP/hSparseRecTrueEff"), kstar.M(), kstar.Pt(), cosRec, centrality, 1.0);
            }
            if (isTestSample) {
              histos.fill(HIST("closureKstarRecoEP/hSparseRecTrue"), kstar.M(), kstar.Pt(), cosRec, centrality, weight);
              histos.fill(HIST("closureKstarRecoEP/hCosThetaStarGenVsRec"), cosGen, cosRec);
            }
          }
          if (!isTestSample) {
            continue;
          }
          histos.fill(HIST("closureKstarRecoEP/hSparseSE"), kstar.M(), kstar.Pt(), cosRec, centrality, weight);

          // rotated background, same rotation as in processSE (weight 1)
          if (closure.fillRotBkg && nBkgRotations.value > 1) {
            const double angleStart = confMinRot.value;
            const double angleStep = (confMaxRot.value - angleStart) / (nBkgRotations.value - 1.0);
            for (int irot = 0; irot < nBkgRotations.value; irot++) {
              const double rotAngle = angleStart + irot * angleStep;
              kaonRot = ROOT::Math::PxPyPzMVector(track1.px() * std::cos(rotAngle) - track1.py() * std::sin(rotAngle), track1.px() * std::sin(rotAngle) + track1.py() * std::cos(rotAngle), track1.pz(), MassKa);
              kstarRot = kaonRot + pion;
              if (std::abs(kstarRot.Rapidity()) > confRapidity) {
                continue;
              }
              if (isPhiReflection(kaonRot, pion)) {
                continue;
              }
              const double cosRot = getCosThetaStar(kstarRot, kaonRot, getSAAxis(kstarRot, psi));
              histos.fill(HIST("closureKstarRecoEP/hSparseRot"), kstarRot.M(), kstarRot.Pt(), cosRot, centrality);
            }
          }
        }
      }

      // ---- generated K*0 / anti-K*0: Acc x Eff denominator (eff sample) and input truth (test sample) ----
      for (const auto& mcParticle : genParticles) {
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (std::abs(mcParticle.pdgCode()) != o2::constants::physics::kK0Star892) {
          continue;
        }
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        if (daughters.size() != NumberOfDaughters) {
          continue;
        }
        bool hasKaon = false;
        bool hasPion = false;
        for (const auto& dau : daughters) {
          if (!dau.isPhysicalPrimary()) {
            continue;
          }
          const bool inAcceptance = !closure.genAcceptanceCut || (dau.pt() > cfgCutPT && std::abs(dau.eta()) < cfgCutEta);
          if (std::abs(dau.pdgCode()) == PDG_t::kKPlus) {
            hasKaon = hasKaon || inAcceptance;
            kaonGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          } else if (std::abs(dau.pdgCode()) == PDG_t::kPiPlus) {
            hasPion = hasPion || inAcceptance;
            pionGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassPi);
          }
        }
        if (!hasKaon || !hasPion) {
          continue;
        }
        kstarGen = kaonGen + pionGen;
        if (std::abs(kstarGen.Rapidity()) > confRapidity) {
          continue;
        }
        const double cosGen = getCosThetaStar(kstarGen, kaonGen, getSAAxis(kstarGen, psiTrue));
        // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
        if (isEffSample) {
          histos.fill(HIST("closureKstarRecoEP/hSparseGenEff"), kstarGen.M(), kstarGen.Pt(), cosGen, centrality, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closureKstarRecoEP/hSparseGen"), kstarGen.M(), kstarGen.Pt(), cosGen, centrality, weight);
        }
      }
    } // rec collision loop
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureRecoEP, "Process MC closure of cos(theta*) for K*, reconstructed Psi_FT0C: SE, true rec, gen", false);

  // ---------------- K*0: mixed event of the test sample (reconstructed Psi_FT0C) ----------------
  void processMCClosureMERecoEP(CollisionMCRecMixTableRecoEP const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    ROOT::Math::PxPyPzMVector kaon, pion, kstar;
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeMCClosureRecoEP binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    SameKindPair<CollisionMCRecMixTableRecoEP, FilTrackMCRecTable, BinningTypeMCClosureRecoEP> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!selectionClosurePairMERecoEP(collision1, collision2, nRecPerMc)) {
        continue;
      }
      const auto centrality = collision1.centFT0C();
      const double psi = getSAEventAngle(collision1.psiFT0C());
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        if (!selectionClosurePrimary(track1) || !selectionClosurePrimary(track2)) {
          continue;
        }
        if (!isTOFOnly && (!strategySelectionPID(track1, 0, strategyPID) || !strategySelectionPID(track2, 1, strategyPID))) {
          continue;
        }
        if (isTOFOnly && (!selectionPID2(track1, 0) || !selectionPID2(track2, 1))) {
          continue;
        }
        kaon = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
        pion = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
        kstar = kaon + pion;
        if (std::abs(kstar.Rapidity()) > confRapidity) {
          continue;
        }
        if (isPhiReflection(kaon, pion)) {
          continue;
        }
        const double cosRec = getCosThetaStar(kstar, kaon, getSAAxis(kstar, psi));
        histos.fill(HIST("closureKstarRecoEP/hSparseME"), kstar.M(), kstar.Pt(), cosRec, centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureMERecoEP, "Process MC closure of cos(theta*) for K*, reconstructed Psi_FT0C: mixed event", false);

  // ---------------- phi(1020): same event, true reco, generated (reconstructed Psi_FT0C) ----------------
  void processMCClosureRecoEPPhi(CollisionMCTrueTable::iterator const& trueCollision, CollisionMCRecTableClosure const& recCollisions, TrackMCTrueTable const& genParticles, FilTrackMCRecTable const& recTracks)
  {
    ROOT::Math::PxPyPzMVector kaonPlus, kaonMinus, phi, kaonPlusRot, phiRot, kaonPlusGen, kaonMinusGen, phiGen;
    static constexpr std::size_t NumberOfDaughters = 2;

    histos.fill(HIST("closurePhiRecoEP/hMC"), 0.5);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("closurePhiRecoEP/hMC"), 1.5);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("closurePhiRecoEP/hMC"), 2.5);
      return;
    }
    const bool isEffSample = isClosureEffSample(trueCollision.globalIndex());
    const bool isTestSample = isClosureTestSample(trueCollision.globalIndex());
    // true axis: defines the injected weight and the Gen/GenEff truth, so the test is not tautological against the smeared axis
    const double psiTrue = getSAEventAngle(trueCollision.eventPlaneAngle());
    for (const auto& recCollision : recCollisions) {
      if (!selectionEventMCRecoEP(recCollision)) {
        continue;
      }
      if (isEffSample) {
        histos.fill(HIST("closurePhiRecoEP/hMC"), 3.5);
      }
      if (isTestSample) {
        histos.fill(HIST("closurePhiRecoEP/hMC"), 4.5);
      }
      const auto centrality = recCollision.centFT0C();
      const double recoPsiFT0C = recCollision.psiFT0C();
      const double psi = getSAEventAngle(recoPsiFT0C); // measured axis: defines cosRec (SE/ME/RecTrue/RecTrueEff/Rot) only
      {
        const double psiTPC = recCollision.psiTPC();
        const double psiFT0A = recCollision.psiFT0A();
        const double qFT0C = recCollision.qFT0C();
        const double qTPC = recCollision.qTPC();
        const double qFT0A = recCollision.qFT0A();
        histos.fill(HIST("closurePhiRecoEP/ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        // plain (unweighted) sub-event resolution, as needed by the EP-method correction formula; same definition as the data ResFT0CTPC etc.
        histos.fill(HIST("closurePhiRecoEP/ResFT0CTPC"), centrality, std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0CFT0A"), centrality, std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
        // closure-only: Psi_FT0C vs the MC true event plane, to cross-check R from sub-events against the true R
        histos.fill(HIST("closurePhiRecoEP/ResTrueFT0C"), centrality, std::cos(2.0 * (recoPsiFT0C - trueCollision.eventPlaneAngle())));
      }
      auto recTracksThisColl = recTracks.sliceBy(perCollision, recCollision.globalIndex());

      // ---- reconstructed pairs ----
      for (const auto& track1 : recTracksThisColl) {
        if (!(track1.signed1Pt() > cfgCutCharge.value)) { // K+
          continue;
        }
        if (!selectionTrackPhi(track1) || !selectionPIDPhi(track1) || !selectionClosurePrimary(track1)) {
          continue;
        }
        for (const auto& track2 : recTracksThisColl) {
          if (!(track2.signed1Pt() < cfgCutCharge.value)) { // K-
            continue;
          }
          if (!selectionTrackPhi(track2) || !selectionPIDPhi(track2) || !selectionClosurePrimary(track2)) {
            continue;
          }
          if (!selectionPairPhi(track1, track2)) {
            continue;
          }
          if (phiSA.removeFakeTrack && (isFakeKaonPhi(track1) || isFakeKaonPhi(track2))) {
            continue;
          }
          kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
          phi = kaonPlus + kaonMinus;
          if (std::abs(phi.Rapidity()) > phiSA.confRapidity) {
            continue;
          }
          const double cosRec = getCosThetaStar(phi, kaonMinus, getSAAxis(phi, psi));

          double weight = 1.0;
          if (isTruePhiPair(track1, track2)) {
            const auto mcPlus = track1.mcParticle();
            const auto mcMinus = track2.mcParticle();
            kaonPlusGen = ROOT::Math::PxPyPzMVector(mcPlus.px(), mcPlus.py(), mcPlus.pz(), MassKa);
            kaonMinusGen = ROOT::Math::PxPyPzMVector(mcMinus.px(), mcMinus.py(), mcMinus.pz(), MassKa);
            phiGen = kaonPlusGen + kaonMinusGen;
            const double cosGen = getCosThetaStar(phiGen, kaonMinusGen, getSAAxis(phiGen, psiTrue));
            weight = getInjectionWeight(cosGen);
            // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
            if (isEffSample) {
              histos.fill(HIST("closurePhiRecoEP/hSparseRecTrueEff"), phi.M(), phi.Pt(), cosRec, centrality, 1.0);
            }
            if (isTestSample) {
              histos.fill(HIST("closurePhiRecoEP/hSparseRecTrue"), phi.M(), phi.Pt(), cosRec, centrality, weight);
              histos.fill(HIST("closurePhiRecoEP/hCosThetaStarGenVsRec"), cosGen, cosRec);
            }
          }
          if (!isTestSample) {
            continue;
          }
          histos.fill(HIST("closurePhiRecoEP/hSparseSE"), phi.M(), phi.Pt(), cosRec, centrality, weight);

          // rotated background: rotate the K+ (weight 1)
          if (closure.fillRotBkg && nBkgRotations.value > 1) {
            const double angleStart = confMinRot.value;
            const double angleStep = (confMaxRot.value - angleStart) / (nBkgRotations.value - 1.0);
            for (int irot = 0; irot < nBkgRotations.value; irot++) {
              const double rotAngle = angleStart + irot * angleStep;
              kaonPlusRot = ROOT::Math::PxPyPzMVector(track1.px() * std::cos(rotAngle) - track1.py() * std::sin(rotAngle), track1.px() * std::sin(rotAngle) + track1.py() * std::cos(rotAngle), track1.pz(), MassKa);
              phiRot = kaonPlusRot + kaonMinus;
              if (std::abs(phiRot.Rapidity()) > phiSA.confRapidity) {
                continue;
              }
              const double cosRot = getCosThetaStar(phiRot, kaonMinus, getSAAxis(phiRot, psi));
              histos.fill(HIST("closurePhiRecoEP/hSparseRot"), phiRot.M(), phiRot.Pt(), cosRot, centrality);
            }
          }
        }
      }

      // ---- generated phi: Acc x Eff denominator (eff sample) and input truth (test sample) ----
      for (const auto& mcParticle : genParticles) {
        if (std::abs(mcParticle.y()) > phiSA.confRapidity) {
          continue;
        }
        if (mcParticle.pdgCode() != o2::constants::physics::kPhi) {
          continue;
        }
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        if (daughters.size() != NumberOfDaughters) {
          continue;
        }
        bool hasKPlus = false;
        bool hasKMinus = false;
        for (const auto& dau : daughters) {
          if (!dau.isPhysicalPrimary()) {
            continue;
          }
          const bool inAcceptance = !closure.genAcceptanceCut || (dau.pt() > cfgCutPT && std::abs(dau.eta()) < cfgCutEta);
          if (dau.pdgCode() == PDG_t::kKPlus) {
            hasKPlus = hasKPlus || inAcceptance;
            kaonPlusGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          } else if (dau.pdgCode() == PDG_t::kKMinus) {
            hasKMinus = hasKMinus || inAcceptance;
            kaonMinusGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), MassKa);
          }
        }
        if (!hasKPlus || !hasKMinus) {
          continue;
        }
        phiGen = kaonPlusGen + kaonMinusGen;
        const double cosGen = getCosThetaStar(phiGen, kaonMinusGen, getSAAxis(phiGen, psiTrue));
        // Acc x Eff must come from an unweighted (unpolarized) sample, independent of the injected rho00 under test
        if (isEffSample) {
          histos.fill(HIST("closurePhiRecoEP/hSparseGenEff"), phiGen.M(), phiGen.Pt(), cosGen, centrality, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closurePhiRecoEP/hSparseGen"), phiGen.M(), phiGen.Pt(), cosGen, centrality, weight);
        }
      }
    } // rec collision loop
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureRecoEPPhi, "Process MC closure of cos(theta*) for phi(1020), reconstructed Psi_FT0C: SE, true rec, gen", false);

  // ---------------- phi(1020): mixed event of the test sample (reconstructed Psi_FT0C) ----------------
  void processMCClosureMERecoEPPhi(CollisionMCRecMixTableRecoEP const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    ROOT::Math::PxPyPzMVector kaonPlus, kaonMinus, phi;
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeMCClosureRecoEP binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    SameKindPair<CollisionMCRecMixTableRecoEP, FilTrackMCRecTable, BinningTypeMCClosureRecoEP> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!selectionClosurePairMERecoEP(collision1, collision2, nRecPerMc)) {
        continue;
      }
      const auto centrality = collision1.centFT0C();
      const double psi = getSAEventAngle(collision1.psiFT0C());
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrackPhi(track1) || !selectionTrackPhi(track2)) {
          continue;
        }
        if (!selectionClosurePrimary(track1) || !selectionClosurePrimary(track2)) {
          continue;
        }
        if (!selectionPIDPhi(track1) || !selectionPIDPhi(track2)) {
          continue;
        }
        if (!selectionPairPhi(track1, track2)) {
          continue;
        }
        if (phiSA.removeFakeTrack && (isFakeKaonPhi(track1) || isFakeKaonPhi(track2))) {
          continue;
        }
        if (track1.sign() > 0) {
          kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
        } else {
          kaonMinus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonPlus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
        }
        phi = kaonPlus + kaonMinus;
        if (std::abs(phi.Rapidity()) > phiSA.confRapidity) {
          continue;
        }
        const double cosRec = getCosThetaStar(phi, kaonMinus, getSAAxis(phi, psi));
        histos.fill(HIST("closurePhiRecoEP/hSparseME"), phi.M(), phi.Pt(), cosRec, centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureMERecoEPPhi, "Process MC closure of cos(theta*) for phi(1020), reconstructed Psi_FT0C: mixed event", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Kstarpbpb>(cfgc)};
}
