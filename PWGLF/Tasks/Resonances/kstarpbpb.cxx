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

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFSpinAlignVecMesonTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
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
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct Kstarpbpb {
  struct RCTCut : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};

    RCTFlagsChecker rctChecker;
  };

  RCTCut rctCut;

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
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
  Configurable<bool> fillLikeSignSA{"fillLikeSignSA", false, "fill the same-event like-sign SA sparse (needs fillSA)"};
  Configurable<bool> useSP{"useSP", false, "use SP"};
  Configurable<int> cfgSAFrame{"cfgSAFrame", 0, "SA quantization axis: 0 = event plane (FT0C), 1 = production plane, 2 = random event plane (random direction in the transverse plane), 3 = random 3D direction (isotropic, per candidate)"};
  Configurable<int> cfgRndSeed{"cfgRndSeed", 0, "Seed for random event plane (0 = unique seed per job)"};
  Configurable<bool> cfgPhiVeto{"cfgPhiVeto", false, "Reject K pi pairs whose KK mass (pion candidate given the kaon mass) is near the phi(1020) mass"};
  Configurable<float> cfgPhiVetoWindow{"cfgPhiVetoWindow", 0.010f, "Half-width of the phi(1020) veto window (GeV/c^2)"};

  struct : ConfigurableGroup {
    std::string prefix = "evtSigLoss";
    Configurable<bool> cutVzGen{"cutVzGen", true, "Cut |vz| < cfgCutVertex on the generated collision"};
    Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Require INEL > 0 for the generated collision"};
    Configurable<bool> isApplyTVX{"isApplyTVX", false, "Require TVX (FT0A and FT0C MC multiplicity > 0) for the generated collision"};
    Configurable<bool> requireSingleReco{"requireSingleReco", true, "Accepted MC event: exactly one reconstructed collision passing the selection, as in processMC / processMCPhi (false: at least one, as in phianalysisrun3pbpb.cxx)"};
    ConfigurableAxis axisImpactPar{"axisImpactPar", {200, 0.0, 20.0}, "Impact parameter (fm)"};
    ConfigurableAxis axisMultMC{"axisMultMC", {500, 0.0, 2000.0}, "Generated charged particles in |eta| < 0.5"};
    ConfigurableAxis axisCosThetaStar{"axisCosThetaStar", {10, -1.0, 1.0}, "cos(#vartheta*) of the generated decay (axis from cfgSAFrame; event plane frame: true reaction plane)"};
  } evtSigLoss;

  struct : ConfigurableGroup {
    std::string prefix = "phiSA";
    Configurable<bool> removeFakeTrack{"removeFakeTrack", true, "Remove fake track from momentum difference"};
    Configurable<float> confFakeKaonCut{"confFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
    Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "Global track (with DCA) + PV contributor + ITS clusters"};
    Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "Global track w/o DCA + PV contributor + manual DCAxy/DCAz + ITS clusters"};
    Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for kaon tracks (manual DCA cut)"};
    Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for kaon tracks (manual DCA cut)"};
    Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
    Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
    Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the combined TPC-TOF Nsigma cut"};
    Configurable<bool> isNoTOF{"isNoTOF", false, "TPC-only PID for all tracks"};
    Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
    Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
    Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
    Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
    Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
    Configurable<bool> fillRotation{"fillRotation", false, "Fill the same-event rotated K+K- background (K+ rotated, with confMinRot, confMaxRot and nBkgRotations as for K*)"};
    Configurable<bool> fillLikeSign{"fillLikeSign", false, "Fill the same-event like-sign (K+K+ and K-K-) background"};
    ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
    ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, -1.0, 1.}, "cos(#vartheta)"};
    ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
    ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {8, 0, 0.8}, "Rapidity"};
    ConfigurableAxis configThnAxisSA{"configThnAxisSA", {200, -1, 1}, "SA"};
  } phiSA;

  struct : ConfigurableGroup {
    std::string prefix = "closure";
    Configurable<bool> injectPolarization{"injectPolarization", false, "Reweight the true K*/phi with W(cos theta*) so that the sample has rho00 = rho00Inj"};
    Configurable<float> rho00Inj{"rho00Inj", 1.0f / 3.0f, "Injected rho00, in [0,1] (1/3 = unpolarized)"};
    Configurable<bool> splitSample{"splitSample", true, "Even MC events build Acc x Eff, odd MC events are the test sample (false: same events for both)"};
    Configurable<bool> fillRotBkg{"fillRotBkg", false, "Fill the rotated-pair background (weight 1), cross-check of the ME background"};
    Configurable<bool> genAcceptanceCut{"genAcceptanceCut", false, "Daughter pT/eta cut on the generated K*/phi (false: denominator is the full |y| < cut sample, so the correction is the full Acc x Eff in cos theta*)"};
    Configurable<bool> requirePrimaryTracks{"requirePrimaryTracks", false, "SE, ME and rotated pairs only from tracks matched to physical primaries (as in phianalysisrun3pbpb / kstar892LightIon); false: all tracks, as in data"};
    Configurable<bool> meRequireSingleReco{"meRequireSingleReco", true, "ME: both events must come from an MC collision with exactly one reco collision, as required in the SE (processMCClosure / processMCClosurePhi)"};
    ConfigurableAxis axisCosThetaStar{"axisCosThetaStar", {10, -1.0, 1.0}, "cos(#vartheta*)"};
  } closure;

  struct : ConfigurableGroup {
    std::string prefix = "mcEffEP";
    Configurable<int> psiSource{"psiSource", 1, "Psi of the EP-frame MC efficiency: 0 = fixed Psi = 0 (old: axis -y for K*, in-plane x for phi), 1 = Psi2 from generated charged primaries in [etaMin, etaMax] (thesis), 2 = MC true event plane (mcCollision.eventPlaneAngle())"};
    Configurable<float> etaMin{"etaMin", -3.3f, "eta min of the generated-track Q2 vector (FT0C acceptance)"};
    Configurable<float> etaMax{"etaMax", -2.1f, "eta max of the generated-track Q2 vector (FT0C acceptance)"};
    Configurable<bool> fillPhiVMinusPsi{"fillPhiVMinusPsi", false, "Fill gen/rec Acc x Eff maps vs pT, cos theta* and phi_V - Psi (thesis Fig. B.69 / Table B.1 check)"};
    ConfigurableAxis axisPhiVMinusPsi{"axisPhiVMinusPsi", {6, 0.0, o2::constants::math::PI}, "#varphi_{V} - #Psi"};
  } mcEffEP;

  struct : ConfigurableGroup {
    std::string prefix = "moment";
    Configurable<bool> fillCosThetaStar{"fillCosThetaStar", true, "Also fill mass, pT, centrality, sample, cos(theta*) sparses next to the moments (SE/ME ratio fit, Method 1)"};
    ConfigurableAxis axisCosThetaStar{"axisCosThetaStar", {20, -1.0, 1.0}, "cos(#vartheta*), full range -1 to 1"};
    Configurable<int> replNDraws{"replNDraws", 5, "Replacement ME: number of donor candidates drawn per target candidate"};
    Configurable<int> replCellCapacity{"replCellCapacity", 10, "Replacement ME: donor reservoir size per matching cell (filled by previous events of the same vz and centrality class)"};
    Configurable<float> replPtWidth{"replPtWidth", 0.2f, "Replacement ME: parent pT matching cell width (GeV/c)"};
    Configurable<float> replEtaWidth{"replEtaWidth", 0.2f, "Replacement ME: parent eta matching cell width"};
    Configurable<float> replPhiWidth{"replPhiWidth", static_cast<float>(o2::constants::math::TwoPI / 18.), "Replacement ME: parent phi matching cell width (rad)"};
    ConfigurableAxis replAxisVertex{"replAxisVertex", {5, -10.0, 10.0}, "Replacement ME: vz classes of the donor pools (centrality classes: axisCentrality; never binned in Psi)"};
    Configurable<bool> mcClosure{"mcClosure", false, "processMCClosureRecoEP / processMCClosureRecoEPPhi also fill the moments: test-sample pairs with the injected weight (sample 0) and efficiency-sample true pairs as the MC reference (sample 2)"};
    Configurable<float> mcRefRho00{"mcRefRho00", 1.0f / 3.0f, "MC reference (sample 2): rho00 put in by weight relative to the MC reaction plane (1/3 = unpolarized reference)"};
    Configurable<int> nSubsamples{"nSubsamples", 30, "Number of event subsamples (resampling units of the offline bootstrap); the uncertainty itself is known to about 1/sqrt(2 (n - 1))"};
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisCentrality{"axisCentrality", {8, 0.0, 80.0}, "FT0C centrality (%), AO2D and derived input"};
  } moment;

  struct : ConfigurableGroup {
    std::string prefix = "candidate";
    Configurable<float> kstarMassMin{"kstarMassMin", 0.70f, "Minimum K pi mass (GeV/c^2)"};
    Configurable<float> kstarMassMax{"kstarMassMax", 1.10f, "Maximum K pi mass (GeV/c^2)"};
    Configurable<float> phiMassMin{"phiMassMin", 0.99f, "Minimum K K mass (GeV/c^2)"};
    Configurable<float> phiMassMax{"phiMassMax", 1.07f, "Maximum K K mass (GeV/c^2)"};
    ConfigurableAxis axisMassKstar{"axisMassKstar", {160, 0.7, 1.1}, "#it{M}_{K#pi} (GeV/#it{c}^{2})"};
    ConfigurableAxis axisMassPhi{"axisMassPhi", {80, 0.99, 1.07}, "#it{M}_{KK} (GeV/#it{c}^{2})"};
    Configurable<float> nSigmaTPC{"nSigmaTPC", 3.0f, "Derived input: maximum |n-sigma TPC| of the daughters (assigned species)"};
    Configurable<float> nSigmaTOF{"nSigmaTOF", 3.0f, "Derived input: maximum |n-sigma TOF| of the daughters with TOF (assigned species)"};
  } candCfg;

  struct : ConfigurableGroup {
    std::string prefix = "v2derived";
    ConfigurableAxis axisPt{"axisPt", {100, 0.0, 10.0}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisCentrality{"axisCentrality", {16, 0.0, 80.0}, "FT0C centrality (%): analysis bins and mixing classes (5% in the paper)"};
    ConfigurableAxis axisUQ{"axisUQ", {400, -16.0, 16.0}, "u2 Q2* and mixed-event background term"};
    ConfigurableAxis axisResSP{"axisResSP", {6000, -30.0, 30.0}, "Q_a Q_b* of the SP resolution"};
    ConfigurableAxis mixAxisVertex{"mixAxisVertex", {10, -10.0, 10.0}, "vz classes of the mixing pools"};
    Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of previous events of the same vz and centrality class each event is mixed with"};
  } v2Derived;

  enum MomentSample : int {
    kMomSEUnlike = 0,
    kMomMERepl = 1,
    kMomMCRef = 2,
    kNMomSamples = 3
  };

  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<bool> additionalEvSel1{"additionalEvSel1", true, "Additional evsel1"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<bool> additionalEvSel4{"additionalEvSel4", true, "Additional evsel4"};

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
  using CollisionMCRecTableClosure = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::EPCalibrationTables>>;
  using CollisionMCRecMixTableRecoEP = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs, aod::EvSels, aod::EPCalibrationTables>>;
  using BinningTypeMCClosureRecoEP = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;

  SliceCache cache;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  enum SAFrame : int {
    kEventPlane = 0,
    kProductionPlane = 1,
    kRandomEventPlane = 2,
    kRandom3D = 3
  };
  TRandom3 rndGen;

  void init(o2::framework::InitContext&)
  {
    rctCut.rctChecker.init(
      rctCut.cfgEvtRCTFlagCheckerLabel,
      rctCut.cfgEvtRCTFlagCheckerZDCCheck,
      rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    if (cfgSAFrame.value < kEventPlane || cfgSAFrame.value > kRandom3D) {
      LOGF(fatal, "cfgSAFrame = %d not supported (0 = event plane, 1 = production plane, 2 = random event plane, 3 = random 3D direction)", cfgSAFrame.value);
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
        if (isRandomFrame()) {
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
      histos.add("hMCEffPsi", "Psi of the MC efficiency;centrality (%);#Psi", kTH2F, {configThnAxisCentrality, {36, -o2::constants::math::PI, o2::constants::math::PI}});
      if (mcEffEP.fillPhiVMinusPsi) {
        const AxisSpec phiVMinusPsiAxis{mcEffEP.axisPhiVMinusPsi, "#varphi_{V} - #Psi"};
        histos.add("hSparseKstarMCGenCosThetaStarPhiVPsi", "Generated K*: p_{T}, cos#theta*, #varphi_{V} - #Psi, centrality", HistType::kTHnSparseD, {configThnAxisPt, configThnAxisV2, phiVMinusPsiAxis, configThnAxisCentrality});
        histos.add("hSparseKstarMCRecCosThetaStarPhiVPsi", "Reconstructed K*: p_{T}, cos#theta*, #varphi_{V} - #Psi, centrality", HistType::kTHnSparseD, {configThnAxisPt, configThnAxisV2, phiVMinusPsiAxis, configThnAxisCentrality});
      }
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
      histos.add("QAbefore/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});

      histos.add("QAbefore/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAbefore/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAafter/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAafter/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAafter/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});

      histos.add("QAbefore/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAbefore/TOF_TPC_Mapka_allpi", "TOF + TPC Combined PID for pion;#sigma_{TOF}^{pion};#sigma_{TPC}^{pion}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allpi", "TOF NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAbefore/TPC_Nsigma_allpi", "TPC NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
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

    if (doprocessSEPhi || doprocessMEPhi || doprocessMCPhi) {
      const AxisSpec thnAxisInvMassPhi{phiSA.configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPtPhi{phiSA.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec thnAxisCosThetaStarPhi{phiSA.configThnAxisCosThetaStar, "cos(#vartheta_{OP})"};
      const AxisSpec thnAxisCentralityPhi{phiSA.configThnAxisCentrality, "Centrality (%)"};
      const AxisSpec thnAxisRapidityPhi{phiSA.configThnAxisRapidity, "Rapidity"};
      const AxisSpec thnAxisSAPhi{phiSA.configThnAxisSA, "SA"};
      if (isRandomFrame()) {
        histos.add("phi/hPsiRandom", "Random event plane angle", kTH2F, {centAxis, phiAxis});
      }
      if (doprocessSEPhi) {
        histos.add("phi/ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
        histos.add("phi/ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
        histos.add("phi/ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
        histos.add("phi/hSparseV2SameEventSA", "hSparseV2SameEventSA", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisSAPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        histos.add("phi/hSparseV2SameEventCosThetaStar", "hSparseV2SameEventCosThetaStar", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisCosThetaStarPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        if (phiSA.fillRotation) {
          histos.add("phi/hSparseV2RotationSA", "hSparseV2RotationSA", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisSAPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
          histos.add("phi/hSparseV2RotationCosThetaStar", "hSparseV2RotationCosThetaStar", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisCosThetaStarPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
        }
        if (phiSA.fillLikeSign) {
          histos.add("phi/hSparseV2SameEventLikeSA", "hSparseV2SameEventLikeSA", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisSAPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
          histos.add("phi/hSparseV2SameEventLikeCosThetaStar", "hSparseV2SameEventLikeCosThetaStar", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisPtPhi, thnAxisCosThetaStarPhi, thnAxisRapidityPhi, thnAxisCentralityPhi});
          histos.add("phi/hMassSameEventLikePP", "Same-event like-sign (++) mass", kTH2F, {thnAxisInvMassPhi, thnAxisCentralityPhi});
          histos.add("phi/hMassSameEventLikeNN", "Same-event like-sign (--) mass", kTH2F, {thnAxisInvMassPhi, thnAxisCentralityPhi});
        }
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
        histos.add("phi/hMCEffPsi", "Psi of the MC efficiency;centrality (%);#Psi", kTH2F, {thnAxisCentralityPhi, {36, -o2::constants::math::PI, o2::constants::math::PI}});
        if (mcEffEP.fillPhiVMinusPsi) {
          const AxisSpec phiVMinusPsiAxis{mcEffEP.axisPhiVMinusPsi, "#varphi_{V} - #Psi"};
          histos.add("phi/hSparseV2MCGenCosThetaStarPhiVPsi", "Generated #phi: p_{T}, cos#theta*, #varphi_{V} - #Psi, centrality", HistType::kTHnSparseF, {thnAxisPtPhi, thnAxisCosThetaStarPhi, phiVMinusPsiAxis, thnAxisCentralityPhi});
          histos.add("phi/hSparseV2MCRecCosThetaStarPhiVPsi", "Reconstructed #phi: p_{T}, cos#theta*, #varphi_{V} - #Psi, centrality", HistType::kTHnSparseF, {thnAxisPtPhi, thnAxisCosThetaStarPhi, phiVMinusPsiAxis, thnAxisCentralityPhi});
        }
      }
    }

    if (doprocessV2DerivedSE || doprocessV2DerivedME) {
      if (v2Derived.nMixedEvents.value < 1) {
        LOGF(fatal, "v2derived.nMixedEvents = %d, must be >= 1", v2Derived.nMixedEvents.value);
      }
      const AxisSpec v2MassKstarAxis{candCfg.axisMassKstar, "#it{M}_{K#pi} (GeV/#it{c}^{2})"};
      const AxisSpec v2MassPhiAxis{candCfg.axisMassPhi, "#it{M}_{KK} (GeV/#it{c}^{2})"};
      const AxisSpec v2PtAxis{v2Derived.axisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec v2CentAxis{v2Derived.axisCentrality, "FT0C centrality (%)"};
      if (doprocessV2DerivedSE) {
        const AxisSpec v2UQAxis{v2Derived.axisUQ, "#it{u}_{2}#it{Q}_{2}^{FT0C*}"};
        const AxisSpec v2ResAxis{v2Derived.axisResSP, "#it{Q}_{a}#it{Q}_{b}^{*}"};
        histos.add("v2derived/ResFT0CTPCSP", "SP resolution: Q^{FT0C} Q^{TPC*}", kTH2F, {v2CentAxis, v2ResAxis});
        histos.add("v2derived/ResFT0CFT0ASP", "SP resolution: Q^{FT0C} Q^{FT0A*}", kTH2F, {v2CentAxis, v2ResAxis});
        histos.add("v2derived/ResFT0ATPCSP", "SP resolution: Q^{FT0A} Q^{TPC*}", kTH2F, {v2CentAxis, v2ResAxis});
        histos.add("v2derived/kstar/hSEUnlike", "K* same event, unlike sign: N^{tot}(M) and u_{2}Q_{2}^{*}", HistType::kTHnSparseF, {v2MassKstarAxis, v2PtAxis, v2UQAxis, v2CentAxis});
        histos.add("v2derived/kstar/hSELikePP", "K^{+}#pi^{+} same event", kTH3F, {v2MassKstarAxis, v2PtAxis, v2CentAxis});
        histos.add("v2derived/kstar/hSELikeNN", "K^{-}#pi^{-} same event", kTH3F, {v2MassKstarAxis, v2PtAxis, v2CentAxis});
        histos.add("v2derived/phi/hSEUnlike", "#phi same event, unlike sign: N^{tot}(M) and u_{2}Q_{2}^{*}", HistType::kTHnSparseF, {v2MassPhiAxis, v2PtAxis, v2UQAxis, v2CentAxis});
        histos.add("v2derived/phi/hSELikePP", "K^{+}K^{+} same event", kTH3F, {v2MassPhiAxis, v2PtAxis, v2CentAxis});
        histos.add("v2derived/phi/hSELikeNN", "K^{-}K^{-} same event", kTH3F, {v2MassPhiAxis, v2PtAxis, v2CentAxis});
      }
      if (doprocessV2DerivedME) {
        const AxisSpec v2BkgAxis{v2Derived.axisUQ, "#sum_{i} #it{u}_{2}^{(i)}#it{Q}_{2}^{(i)*} cos 2(#varphi_{i} - #varphi)"};
        histos.add("v2derived/kstar/hMEUnlike", "K* mixed event, unlike sign: N^{mix}(M) and Eq. 5 numerator", HistType::kTHnSparseF, {v2MassKstarAxis, v2PtAxis, v2BkgAxis, v2CentAxis});
        histos.add("v2derived/kstar/hMELikePP", "K^{+}#pi^{+} mixed event", kTH3F, {v2MassKstarAxis, v2PtAxis, v2CentAxis});
        histos.add("v2derived/kstar/hMELikeNN", "K^{-}#pi^{-} mixed event", kTH3F, {v2MassKstarAxis, v2PtAxis, v2CentAxis});
        histos.add("v2derived/phi/hMEUnlike", "#phi mixed event, unlike sign: N^{mix}(M) and Eq. 5 numerator", HistType::kTHnSparseF, {v2MassPhiAxis, v2PtAxis, v2BkgAxis, v2CentAxis});
        histos.add("v2derived/phi/hMELikePP", "K^{+}K^{+} mixed event", kTH3F, {v2MassPhiAxis, v2PtAxis, v2CentAxis});
        histos.add("v2derived/phi/hMELikeNN", "K^{-}K^{-} mixed event", kTH3F, {v2MassPhiAxis, v2PtAxis, v2CentAxis});
        histos.add("v2derived/hMixedEventPairs", "Mixed event pairs;FT0C centrality (%);Pairs", kTH1F, {v2CentAxis});
      }
    }

    const bool replKstar = doprocessMomentDerivedMERepl || doprocessMCClosureReplacementME;
    const bool replPhi = doprocessMomentDerivedMERepl || doprocessMCClosureReplacementMEPhi;
    const bool mcRefKstar = moment.mcClosure && doprocessMCClosureRecoEP;
    const bool mcRefPhi = moment.mcClosure && doprocessMCClosureRecoEPPhi;
    const bool momentKstar = doprocessMomentDerivedSE || replKstar || mcRefKstar;
    const bool momentPhi = doprocessMomentDerivedSE || replPhi || mcRefPhi;
    if (momentKstar || momentPhi) {
      if (moment.nSubsamples.value < 2) {
        LOGF(fatal, "moment.nSubsamples = %d, must be >= 2 for the bootstrap uncertainty", moment.nSubsamples.value);
      }
      const AxisSpec momPtAxis{moment.axisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec momCentAxis{moment.axisCentrality, "Centrality (%)"};
      const AxisSpec momSampleAxis{static_cast<int>(kNMomSamples), -0.5, static_cast<int>(kNMomSamples) - 0.5, "sample (0 SE unlike, 1 ME replacement, 2 MC reference)"};
      const AxisSpec momCosAxis{moment.axisCosThetaStar, "cos(#vartheta*)"};
      const AxisSpec momSubAxis{moment.nSubsamples.value, -0.5, moment.nSubsamples.value - 0.5, "subsample"};
      if (momentKstar) {
        const AxisSpec momMassAxis{candCfg.axisMassKstar, "#it{M}_{K#pi} (GeV/#it{c}^{2})"};
        histos.add("moment/kstar/hX0", "K*: #Sigma 1", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        histos.add("moment/kstar/hX2", "K*: #Sigma cos^{2}#theta*", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        histos.add("moment/kstar/hX4", "K*: #Sigma cos^{4}#theta*", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        if (moment.fillCosThetaStar) {
          histos.add("moment/kstar/hCosThetaStar", "K*: cos#theta* (SE/ME ratio fit)", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momCosAxis});
        }
        if (replKstar) {
          histos.add("moment/kstar/hXD", "K*: #Sigma cos 2(#Psi_{target} - #Psi_{donor}), replacement ME", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        }
        if (doprocessMCClosureReplacementME) {
          histos.add("moment/kstar/hXC", "K*: #Sigma cos 2(#Psi^{EP}_{target} - #Psi^{RP}_{donor}), MC replacement ME", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        }
      }
      if (momentPhi) {
        const AxisSpec momMassAxis{candCfg.axisMassPhi, "#it{M}_{KK} (GeV/#it{c}^{2})"};
        histos.add("moment/phi/hX0", "#phi: #Sigma 1", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        histos.add("moment/phi/hX2", "#phi: #Sigma cos^{2}#theta*", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        histos.add("moment/phi/hX4", "#phi: #Sigma cos^{4}#theta*", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        if (moment.fillCosThetaStar) {
          histos.add("moment/phi/hCosThetaStar", "#phi: cos#theta* (SE/ME ratio fit)", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momCosAxis});
        }
        if (replPhi) {
          histos.add("moment/phi/hXD", "#phi: #Sigma cos 2(#Psi_{target} - #Psi_{donor}), replacement ME", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        }
        if (doprocessMCClosureReplacementMEPhi) {
          histos.add("moment/phi/hXC", "#phi: #Sigma cos 2(#Psi^{EP}_{target} - #Psi^{RP}_{donor}), MC replacement ME", HistType::kTHnSparseD, {momMassAxis, momPtAxis, momCentAxis, momSampleAxis, momSubAxis});
        }
      }
      if (replKstar || replPhi) {
        if (cfgSAFrame.value == kProductionPlane) {
          LOGF(fatal, "The replacement ME projects the donor decay on the target event plane; it needs cfgSAFrame = 0 (event plane), 2 (random event plane) or 3 (random 3D direction)");
        }
        if (moment.replNDraws.value < 1 || moment.replCellCapacity.value < 1 || moment.replPtWidth.value <= 0.f || moment.replEtaWidth.value <= 0.f || moment.replPhiWidth.value <= 0.f) {
          LOGF(fatal, "Replacement ME: replNDraws and replCellCapacity must be >= 1 and the matching cell widths > 0");
        }
        histos.add("moment/hReplCos2DeltaPsi", "Replacement ME: target vs donor plane;centrality (%);cos 2(#Psi_{target} - #Psi_{donor})", kTH2F, {centAxis, {200, -1.0, 1.0}});
        if (replKstar) {
          histos.add("moment/kstar/hReplacement", "K* replacement ME;;target candidates", kTH1F, {{2, -0.5, 1.5}});
          auto hRepl = histos.get<TH1>(HIST("moment/kstar/hReplacement"));
          hRepl->GetXaxis()->SetBinLabel(1, "replaced");
          hRepl->GetXaxis()->SetBinLabel(2, "no donor in cell");
        }
        if (replPhi) {
          histos.add("moment/phi/hReplacement", "#phi replacement ME;;target candidates", kTH1F, {{2, -0.5, 1.5}});
          auto hRepl = histos.get<TH1>(HIST("moment/phi/hReplacement"));
          hRepl->GetXaxis()->SetBinLabel(1, "replaced");
          hRepl->GetXaxis()->SetBinLabel(2, "no donor in cell");
        }
      }
      if (doprocessMomentDerivedSE || mcRefKstar || mcRefPhi) {
        histos.add("moment/ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
        histos.add("moment/ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
        histos.add("moment/ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
      }
      if (mcRefKstar || mcRefPhi) {
        if (moment.mcRefRho00.value < 0.f || moment.mcRefRho00.value > 1.f) {
          LOGF(fatal, "moment.mcRefRho00 = %f is outside [0,1]", moment.mcRefRho00.value);
        }
        histos.add("moment/hResMCTrueRef", "MC reference events;centrality (%);cos 2(#Psi_{EP} - #Psi_{RP})", kTH2F, {centAxis, resAxis});
        histos.add("moment/hResMCTrueTest", "MC test events;centrality (%);cos 2(#Psi_{EP} - #Psi_{RP})", kTH2F, {centAxis, resAxis});
      }
    }

    if (doprocessEvtLossSigLossMC) {
      const AxisSpec impactParAxis{evtSigLoss.axisImpactPar, "Impact parameter (fm)"};
      const AxisSpec multMCAxis{evtSigLoss.axisMultMC, "N_{ch}^{gen} (|#eta| < 0.5)"};
      const AxisSpec centLossAxis{100, 0.0, 100.0, "Centrality FT0C (%)"};
      const AxisSpec ptKstarAxis{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec ptPhiAxis{phiSA.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec cosLossAxis{evtSigLoss.axisCosThetaStar, "cos(#vartheta*)"};
      histos.add("evtSigLoss/MCEventHist", "MC event statistics", kTH1F, {{4, 0.5, 4.5}});
      auto hstat = histos.get<TH1>(HIST("evtSigLoss/MCEventHist"));
      hstat->GetXaxis()->SetBinLabel(1, "All MC events");
      hstat->GetXaxis()->SetBinLabel(2, "MC events with selected reco event");
      hstat->GetXaxis()->SetBinLabel(3, "MC events with no reco event");
      hstat->GetXaxis()->SetBinLabel(4, "MC events with >1 reco event");
      histos.add("evtSigLoss/hImpactParameterGen", "Impact parameter, all generated events", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParameterRec", "Impact parameter, generated events with selected reco event", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParameterGenNoReco", "Impact parameter, generated events with no reco event", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParVsCentRec", "Impact parameter vs centrality of the selected reco event", kTH2F, {centLossAxis, impactParAxis});
      histos.add("evtSigLoss/hMultEta05Gen", "N_{ch}^{gen}, all generated events", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05Rec", "N_{ch}^{gen}, generated events with selected reco event", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05GenNoReco", "N_{ch}^{gen}, generated events with no reco event", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05VsCentRec", "N_{ch}^{gen} vs centrality of the selected reco event", kTH2F, {centLossAxis, multMCAxis});
      histos.add("evtSigLoss/hKstarGenBeforeEvtSel", "Generated K*0, all events", kTH3F, {ptKstarAxis, cosLossAxis, impactParAxis});
      histos.add("evtSigLoss/hKstarGenAfterEvtSel", "Generated K*0, events with selected reco event", kTH3F, {ptKstarAxis, cosLossAxis, impactParAxis});
      histos.add("evtSigLoss/hKstarGenVsMultBeforeEvtSel", "Generated K*0, all events", kTH3F, {ptKstarAxis, cosLossAxis, multMCAxis});
      histos.add("evtSigLoss/hKstarGenVsMultAfterEvtSel", "Generated K*0, events with selected reco event", kTH3F, {ptKstarAxis, cosLossAxis, multMCAxis});
      histos.add("evtSigLoss/hKstarGenAfterEvtSelVsCent", "Generated K*0, events with selected reco event", kTH3F, {ptKstarAxis, cosLossAxis, centLossAxis});
      histos.add("evtSigLoss/hPhiGenBeforeEvtSel", "Generated phi, all events", kTH3F, {ptPhiAxis, cosLossAxis, impactParAxis});
      histos.add("evtSigLoss/hPhiGenAfterEvtSel", "Generated phi, events with selected reco event", kTH3F, {ptPhiAxis, cosLossAxis, impactParAxis});
      histos.add("evtSigLoss/hPhiGenVsMultBeforeEvtSel", "Generated phi, all events", kTH3F, {ptPhiAxis, cosLossAxis, multMCAxis});
      histos.add("evtSigLoss/hPhiGenVsMultAfterEvtSel", "Generated phi, events with selected reco event", kTH3F, {ptPhiAxis, cosLossAxis, multMCAxis});
      histos.add("evtSigLoss/hPhiGenAfterEvtSelVsCent", "Generated phi, events with selected reco event", kTH3F, {ptPhiAxis, cosLossAxis, centLossAxis});
    }

    if (doprocessMCClosureRecoEP || doprocessMCClosureMERecoEP || doprocessMCClosureRecoEPPhi || doprocessMCClosureMERecoEPPhi || doprocessMCClosureReplacementME || doprocessMCClosureReplacementMEPhi) {
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
          histos.add("closureKstarRecoEP/hMC", "MC event statistics", kTH1F, {{5, 0.0f, 5.0f}});
          histos.add("closureKstarRecoEP/hSparseGenEff", "Generated K*, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstarRecoEP/hSparseRecTrueEff", "Reconstructed true K*, efficiency sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstarRecoEP/hSparseGen", "Generated K* (input), test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstarRecoEP/hSparseRecTrue", "Reconstructed true K*, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          histos.add("closureKstarRecoEP/hSparseSE", "Same-event unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          if (closure.fillRotBkg) {
            histos.add("closureKstarRecoEP/hSparseRot", "Rotated unlike-sign K#pi pairs, test sample", HistType::kTHnSparseD, {massAxis, ptAxis, cosAxisClosure, centAxisClosureKstar}, true);
          }
          histos.add("closureKstarRecoEP/hCosThetaStarGenVsRec", "cos(#vartheta*) gen vs rec, true K*", kTH2F, {cosFineAxis, cosFineAxis});
          histos.add("closureKstarRecoEP/ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
          histos.add("closureKstarRecoEP/ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
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
    while (result >= o2::constants::math::PI) {
      result -= o2::constants::math::PI;
    }
    return result;
  }

  bool isRandomFrame() const
  {
    return cfgSAFrame.value == kRandomEventPlane || cfgSAFrame.value == kRandom3D;
  }

  ROOT::Math::XYZVector getRandom3DAxis()
  {
    const double cosTheta = rndGen.Uniform(-1., 1.);
    const double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
    const double phi = rndGen.Uniform(0., o2::constants::math::TwoPI);
    return {sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
  }

  double getSAEventAngle(double psiEP)
  {
    if (isRandomFrame()) {
      return rndGen.Uniform(-o2::constants::math::PIHalf, o2::constants::math::PIHalf);
    }
    return psiEP;
  }

  double getSAPlaneAngle(const ROOT::Math::PxPyPzMVector& mother, double psiSA)
  {
    if (cfgSAFrame.value == kProductionPlane) {
      return mother.Phi();
    }
    return psiSA;
  }

  bool isPhiReflection(const ROOT::Math::PxPyPzMVector& kaon, const ROOT::Math::PxPyPzMVector& pion)
  {
    if (!cfgPhiVeto) {
      return false;
    }
    ROOT::Math::PxPyPzMVector kaonAsKaon(kaon.Px(), kaon.Py(), kaon.Pz(), MassKa);
    ROOT::Math::PxPyPzMVector pionAsKaon(pion.Px(), pion.Py(), pion.Pz(), MassKa);
    return std::abs((kaonAsKaon + pionAsKaon).M() - o2::constants::physics::MassPhi) < cfgPhiVetoWindow;
  }

  ROOT::Math::XYZVector getEPAxis(double psiSA)
  {
    return {std::sin(psiSA), -std::cos(psiSA), 0.};
  }

  ROOT::Math::XYZVector getSAAxis(const ROOT::Math::PxPyPzMVector& mother, double psiSA)
  {
    if (cfgSAFrame.value == kProductionPlane) {
      return ROOT::Math::XYZVector(0., 0., 1.).Cross(mother.Vect()).Unit();
    }
    if (cfgSAFrame.value == kRandom3D) {
      return getRandom3DAxis();
    }
    return getEPAxis(psiSA);
  }

  int getMomentSubsample(int64_t eventIndex)
  {
    return static_cast<int>(eventIndex % moment.nSubsamples.value);
  }

  template <bool IsPhi>
  void fillMoments(double mass, double pt, double centrality, int sample, int subsample, double cosThetaStar, double weight = 1.0)
  {
    const double x2 = cosThetaStar * cosThetaStar;
    const double x4 = x2 * x2;
    if constexpr (IsPhi) {
      histos.fill(HIST("moment/phi/hX0"), mass, pt, centrality, sample, subsample, weight);
      histos.fill(HIST("moment/phi/hX2"), mass, pt, centrality, sample, subsample, weight * x2);
      histos.fill(HIST("moment/phi/hX4"), mass, pt, centrality, sample, subsample, weight * x4);
      if (moment.fillCosThetaStar) {
        histos.fill(HIST("moment/phi/hCosThetaStar"), mass, pt, centrality, sample, cosThetaStar, weight);
      }
    } else {
      histos.fill(HIST("moment/kstar/hX0"), mass, pt, centrality, sample, subsample, weight);
      histos.fill(HIST("moment/kstar/hX2"), mass, pt, centrality, sample, subsample, weight * x2);
      histos.fill(HIST("moment/kstar/hX4"), mass, pt, centrality, sample, subsample, weight * x4);
      if (moment.fillCosThetaStar) {
        histos.fill(HIST("moment/kstar/hCosThetaStar"), mass, pt, centrality, sample, cosThetaStar, weight);
      }
    }
  }

  template <typename TCollision>
  bool selectionEventRCT(const TCollision& collision)
  {
    return !rctCut.requireRCTFlagChecker || rctCut.rctChecker(collision);
  }

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

  template <typename TCollision>
  bool selectionEventCuts(const TCollision& collision)
  {
    return collision.sel8() && collision.triggereventep() && selectionEventBits(collision);
  }

  template <typename TCollision>
  bool selectionEvent(const TCollision& collision)
  {
    return selectionEventRCT(collision) && selectionEventCuts(collision);
  }

  template <typename TCollision>
  bool selectionEventPairME(const TCollision& collision1, const TCollision& collision2)
  {
    if (!selectionEvent(collision1) || !selectionEvent(collision2)) {
      return false;
    }
    return collision1.bcId() != collision2.bcId();
  }

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

  template <typename T1, typename T2>
  bool selectionPairPhi(const T1& candidate1, const T2& candidate2)
  {
    double pt1 = candidate1.pt(), pt2 = candidate2.pt();
    double pz1 = candidate1.pz(), pz2 = candidate2.pz();
    double p1 = candidate1.p(), p2 = candidate2.p();
    double angle = std::acos(std::clamp((pt1 * pt2 + pz1 * pz2) / (p1 * p2), -1.0, 1.0));
    return !phiSA.isDeepAngle || angle >= phiSA.cfgDeepAngle;
  }

  template <typename T>
  bool isFakeKaonPhi(T const& track)
  {
    return std::abs(track.p() - track.tpcInnerParam()) > phiSA.confFakeKaonCut;
  }

  ROOT::Math::XYZVector getSAAxisPhiMC(const ROOT::Math::PxPyPzMVector& mother, double psiSA)
  {
    if (cfgSAFrame.value == kEventPlane && mcEffEP.psiSource == 0) {
      return {std::cos(2.0 * psiSA), std::sin(2.0 * psiSA), 0.};
    }
    return getSAAxis(mother, psiSA);
  }

  template <typename TParticles>
  double getGenPsi2(TParticles const& particles)
  {
    double qx = 0., qy = 0.;
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary() || particle.eta() < mcEffEP.etaMin || particle.eta() > mcEffEP.etaMax) {
        continue;
      }
      const int absPdg = std::abs(particle.pdgCode());
      if (absPdg != PDG_t::kPiPlus && absPdg != PDG_t::kKPlus && absPdg != PDG_t::kProton && absPdg != PDG_t::kElectron && absPdg != PDG_t::kMuonMinus) {
        continue;
      }
      qx += std::cos(2.0 * particle.phi());
      qy += std::sin(2.0 * particle.phi());
    }
    return 0.5 * std::atan2(qy, qx);
  }

  template <typename TMcCollision, typename TParticles>
  double getMCEffPsi(TMcCollision const& mcCollision, TParticles const& particles)
  {
    if (mcEffEP.psiSource == 1) {
      return getGenPsi2(particles);
    }
    if (mcEffEP.psiSource == 2) {
      return mcCollision.eventPlaneAngle();
    }
    return 0.0;
  }

  static constexpr float HalfPI = o2::constants::math::PI * 0.5f;
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle",
                               {6, -HalfPI, HalfPI},
                               "event plane angle"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
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
    auto psiSA = getSAEventAngle(psiFT0C);
    if (fillSA && isRandomFrame()) {
      histos.fill(HIST("hPsiRandom"), centrality, psiSA);
    }
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
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
        kstarMother = daughter1 + daughter2;
        if (std::abs(kstarMother.Rapidity()) > confRapidity) {
          continue;
        }
        if (isPhiReflection(daughter1, daughter2)) {
          continue;
        }
        auto phiMinusPsi = getPhiInRange(kstarMother.Phi() - psiFT0C);

        if (useSP) {
          v2 = std::cos(2.0 * phiMinusPsi) * qFT0C;
        }
        if (!useSP) {
          v2 = std::cos(2.0 * phiMinusPsi);
        }
        if (additionalQAplots1) {
          histos.fill(HIST("ResTrackSPFT0CTPC"), centrality, qFT0C * qTPC * std::cos(2.0 * (psiFT0C - psiTPC)));
          histos.fill(HIST("ResTrackSPFT0CFT0A"), centrality, qFT0C * qFT0A * std::cos(2.0 * (psiFT0C - psiFT0A)));
          histos.fill(HIST("ResTrackSPFT0ATPC"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        }
        if (!fillSA) {
          histos.fill(HIST("hSparseV2SASameEvent_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
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
            if (isPhiReflection(kaonrot, daughter2)) {
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
  PROCESS_SWITCH(Kstarpbpb, processSE, "K* AO2D: same event (SA or v2 sparses)", true);

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    ROOT::Math::PxPyPzMVector kstarMother, fourVecDauCM, daughter1, daughter2;
    ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVecNorm;
    double v2 = 0.;

    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}};
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
        if (isPhiReflection(daughter1, daughter2)) {
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
  PROCESS_SWITCH(Kstarpbpb, processMixedEvent, "K* AO2D: mixed event (combinatorial background, SA or v2 sparses)", true);

  void processMC(CollisionMCTrueTable::iterator const& TrueCollision, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
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
      auto psiFT0C = getMCEffPsi(TrueCollision, GenParticles);
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
      auto psiSA = getSAEventAngle(psiFT0C);
      histos.fill(HIST("hMCEffPsi"), centrality, psiSA);
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
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
              kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
              pionMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
              kstarMother = kaonPlus + pionMinus;
              if (std::abs(kstarMother.Rapidity()) > confRapidity) {
                continue;
              }
              if (isPhiReflection(kaonPlus, pionMinus)) {
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
                if (mcEffEP.fillPhiVMinusPsi) {
                  histos.fill(HIST("hSparseKstarMCRecCosThetaStarPhiVPsi"), kstarMother.Pt(), cosThetaStarRec, getPhiInRange(kstarMother.Phi() - psiSA), centrality);
                }
              }
            }
          }
        }
      }
      for (const auto& mcParticle : GenParticles) {
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (pdgcheck && std::abs(mcParticle.pdgCode()) != o2::constants::physics::kK0Star892) {
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
            if (mcEffEP.fillPhiVMinusPsi) {
              histos.fill(HIST("hSparseKstarMCGenCosThetaStarPhiVPsi"), kstarMother.Pt(), cosThetaStarGen, getPhiInRange(kstarMother.Phi() - psiSA), centrality);
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMC, "K* MC: Acc x Eff (gen, rec)", false);

  struct SAValues {
    double cosThetaStar;
    double sa;
  };

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
    if (isRandomFrame()) {
      histos.fill(HIST("phi/hPsiRandom"), centrality, psiSA);
    }
    histos.fill(HIST("phi/ResFT0CTPC"), centrality, std::cos(2.0 * (collision.psiFT0C() - collision.psiTPC())));
    histos.fill(HIST("phi/ResFT0CFT0A"), centrality, std::cos(2.0 * (collision.psiFT0C() - collision.psiFT0A())));
    histos.fill(HIST("phi/ResFT0ATPC"), centrality, std::cos(2.0 * (collision.psiTPC() - collision.psiFT0A())));

    ROOT::Math::PxPyPzMVector kaonPlusPhi, kaonMinusPhi, phiMother, kaonPlusRot, phiRot;
    for (const auto& track1 : tracks) {
      if (!(track1.signed1Pt() > cfgCutCharge.value)) {
        continue;
      }
      if (!selectionTrackPhi(track1) || !selectionPIDPhi(track1)) {
        continue;
      }
      for (const auto& track2 : tracks) {
        if (!(track2.signed1Pt() < cfgCutCharge.value)) {
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
        histos.fill(HIST("phi/hSparseV2SameEventSA"), phiMother.M(), phiMother.Pt(), sa, absRapidity, centrality);
        histos.fill(HIST("phi/hSparseV2SameEventCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, absRapidity, centrality);
        if (phiSA.fillRotation && nBkgRotations.value > 1) {
          const double angleStart = confMinRot.value;
          const double angleStep = (confMaxRot.value - angleStart) / (nBkgRotations.value - 1.0);
          for (int irot = 0; irot < nBkgRotations.value; irot++) {
            const double rotAngle = angleStart + irot * angleStep;
            kaonPlusRot = ROOT::Math::PxPyPzMVector(track1.px() * std::cos(rotAngle) - track1.py() * std::sin(rotAngle), track1.px() * std::sin(rotAngle) + track1.py() * std::cos(rotAngle), track1.pz(), MassKa);
            phiRot = kaonPlusRot + kaonMinusPhi;
            auto absRapidityRot = std::abs(phiRot.Rapidity());
            if (absRapidityRot > phiSA.confRapidity) {
              continue;
            }
            auto [cosThetaStarRot, saRot] = getSAValuesPhi(phiRot, kaonMinusPhi, getSAAxis(phiRot, psiSA), psiSA);
            histos.fill(HIST("phi/hSparseV2RotationSA"), phiRot.M(), phiRot.Pt(), saRot, absRapidityRot, centrality);
            histos.fill(HIST("phi/hSparseV2RotationCosThetaStar"), phiRot.M(), phiRot.Pt(), cosThetaStarRot, absRapidityRot, centrality);
          }
        }
      }
    }

    if (phiSA.fillLikeSign) {
      ROOT::Math::PxPyPzMVector kaonLike1, kaonLike2, pairLike;
      for (const auto& track1 : tracks) {
        if (!selectionTrackPhi(track1) || !selectionPIDPhi(track1)) {
          continue;
        }
        for (const auto& track2 : tracks) {
          if (track2.globalIndex() <= track1.globalIndex()) {
            continue;
          }
          if (track1.sign() * track2.sign() <= 0) {
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
          kaonLike1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
          kaonLike2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassKa);
          pairLike = kaonLike1 + kaonLike2;
          auto absRapidityLike = std::abs(pairLike.Rapidity());
          if (absRapidityLike > phiSA.confRapidity) {
            continue;
          }
          auto [cosThetaStarLike, saLike] = getSAValuesPhi(pairLike, kaonLike1, getSAAxis(pairLike, psiSA), psiSA);
          histos.fill(HIST("phi/hSparseV2SameEventLikeSA"), pairLike.M(), pairLike.Pt(), saLike, absRapidityLike, centrality);
          histos.fill(HIST("phi/hSparseV2SameEventLikeCosThetaStar"), pairLike.M(), pairLike.Pt(), cosThetaStarLike, absRapidityLike, centrality);
          if (track1.sign() > 0) {
            histos.fill(HIST("phi/hMassSameEventLikePP"), pairLike.M(), centrality);
          } else {
            histos.fill(HIST("phi/hMassSameEventLikeNN"), pairLike.M(), centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processSEPhi, "phi AO2D: same event (SA sparses)", false);

  void processMEPhi(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}};
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
        histos.fill(HIST("phi/hSparseV2MixedEventSA"), phiMother.M(), phiMother.Pt(), sa, absRapidity, centrality);
        histos.fill(HIST("phi/hSparseV2MixedEventCosThetaStar"), phiMother.M(), phiMother.Pt(), cosThetaStar, absRapidity, centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMEPhi, "phi AO2D: mixed event (combinatorial background)", false);

  void processMCPhi(CollisionMCTrueTable::iterator const& TrueCollision, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
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
      auto psiFT0C = getMCEffPsi(TrueCollision, GenParticles);
      histos.fill(HIST("phi/hMC"), 3);
      if (!RecCollision.sel8()) {
        histos.fill(HIST("phi/hMC"), 4);
        continue;
      }
      if (!selectionEventBits(RecCollision)) {
        continue;
      }
      histos.fill(HIST("phi/hMC"), 5);
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("phi/hMC"), 6);
        continue;
      }
      histos.fill(HIST("phi/hMC"), 7);
      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("phi/CentPercentileMCRecHist"), centrality);
      auto psiSA = getSAEventAngle(psiFT0C);
      histos.fill(HIST("phi/hMCEffPsi"), centrality, psiSA);
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
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
              if (std::abs(phiMother.Rapidity()) > phiSA.confRapidity) {
                continue;
              }
              auto [cosThetaStar, sa] = getSAValuesPhi(phiMother, kaonMinusPhi, getSAAxisPhiMC(phiMother, psiSA), psiSA);
              histos.fill(HIST("phi/hSparseV2MCRecCosThetaStar_effy"), phiMother.M(), phiMother.Pt(), cosThetaStar, std::abs(phiMother.Rapidity()), centrality);
              if (mcEffEP.fillPhiVMinusPsi) {
                histos.fill(HIST("phi/hSparseV2MCRecCosThetaStarPhiVPsi"), phiMother.Pt(), cosThetaStar, getPhiInRange(phiMother.Phi() - psiSA), centrality);
              }
              histos.fill(HIST("phi/hSparseV2MCRecSA"), phiMother.M(), phiMother.Pt(), sa, std::abs(phiMother.Rapidity()), centrality);
            }
          }
        }
      }
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
          if (mcEffEP.fillPhiVMinusPsi) {
            histos.fill(HIST("phi/hSparseV2MCGenCosThetaStarPhiVPsi"), phiMother.Pt(), cosThetaStar, getPhiInRange(phiMother.Phi() - psiSA), centrality);
          }
          histos.fill(HIST("phi/hSparseV2MCGenSA"), phiMother.M(), phiMother.Pt(), sa, std::abs(phiMother.Rapidity()), centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMCPhi, "phi MC: Acc x Eff (gen, rec)", false);

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

    static constexpr std::size_t NumberOfDaughters = 2;
    const double psiTrue = getSAEventAngle(mcCollision.eventPlaneAngle());
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
      const ROOT::Math::PxPyPzMVector mother = dau1 + dau2;
      const double pt = mother.Pt();
      if (isKstar && hasKaon && hasPion) {
        const double cosThetaStar = getCosThetaStar(mother, dau1, getSAAxis(mother, psiTrue));
        histos.fill(HIST("evtSigLoss/hKstarGenBeforeEvtSel"), pt, cosThetaStar, impactPar);
        histos.fill(HIST("evtSigLoss/hKstarGenVsMultBeforeEvtSel"), pt, cosThetaStar, multMC);
        if (selected) {
          histos.fill(HIST("evtSigLoss/hKstarGenAfterEvtSel"), pt, cosThetaStar, impactPar);
          histos.fill(HIST("evtSigLoss/hKstarGenVsMultAfterEvtSel"), pt, cosThetaStar, multMC);
          histos.fill(HIST("evtSigLoss/hKstarGenAfterEvtSelVsCent"), pt, cosThetaStar, centrality);
        }
      } else if (!isKstar && hasKPlus && hasKMinus) {
        const double cosThetaStar = getCosThetaStar(mother, dau2, getSAAxis(mother, psiTrue));
        histos.fill(HIST("evtSigLoss/hPhiGenBeforeEvtSel"), pt, cosThetaStar, impactPar);
        histos.fill(HIST("evtSigLoss/hPhiGenVsMultBeforeEvtSel"), pt, cosThetaStar, multMC);
        if (selected) {
          histos.fill(HIST("evtSigLoss/hPhiGenAfterEvtSel"), pt, cosThetaStar, impactPar);
          histos.fill(HIST("evtSigLoss/hPhiGenVsMultAfterEvtSel"), pt, cosThetaStar, multMC);
          histos.fill(HIST("evtSigLoss/hPhiGenAfterEvtSelVsCent"), pt, cosThetaStar, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processEvtLossSigLossMC, "K* and phi MC: event loss and signal loss", false);

  double getRho00Weight(double cosThetaStar, double rho)
  {
    return 1.5 * ((1.0 - rho) + (3.0 * rho - 1.0) * cosThetaStar * cosThetaStar);
  }

  double getInjectionWeight(double cosThetaStar)
  {
    if (!closure.injectPolarization) {
      return 1.0;
    }
    return getRho00Weight(cosThetaStar, closure.rho00Inj.value);
  }

  double getCosThetaStar(const ROOT::Math::PxPyPzMVector& mother, const ROOT::Math::PxPyPzMVector& dau, const ROOT::Math::XYZVector& axis)
  {
    ROOT::Math::Boost boost{mother.BoostToCM()};
    auto threeVecDau = boost(dau).Vect();
    return axis.Dot(threeVecDau) / std::sqrt(threeVecDau.Mag2()) / std::sqrt(axis.Mag2());
  }

  static constexpr int64_t SplitSampleModulus = 2;
  bool isClosureEffSample(int64_t mcCollisionIndex)
  {
    return !closure.splitSample || (mcCollisionIndex % SplitSampleModulus == 0);
  }
  bool isClosureTestSample(int64_t mcCollisionIndex)
  {
    return !closure.splitSample || (mcCollisionIndex % SplitSampleModulus == 1);
  }

  template <typename TTrack>
  bool selectionClosurePrimary(const TTrack& track)
  {
    return !closure.requirePrimaryTracks || (track.has_mcParticle() && track.mcParticle().isPhysicalPrimary());
  }

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

  template <typename TCollision>
  bool selectionEventMC(const TCollision& collision)
  {
    return collision.sel8() && selectionEventBits(collision) && std::abs(collision.posZ()) < cfgCutVertex && std::abs(collision.centFT0C()) < cfgCutCentrality;
  }

  template <typename TCollision>
  bool selectionEventMCRecoEP(const TCollision& collision)
  {
    return selectionEventMC(collision) && collision.triggereventep();
  }

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
      const double psi = getSAEventAngle(recoPsiFT0C);
      {
        const double psiTPC = recCollision.psiTPC();
        const double psiFT0A = recCollision.psiFT0A();
        const double qFT0C = recCollision.qFT0C();
        const double qTPC = recCollision.qTPC();
        const double qFT0A = recCollision.qFT0A();
        histos.fill(HIST("closureKstarRecoEP/ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0CTPC"), centrality, std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0CFT0A"), centrality, std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closureKstarRecoEP/ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
        histos.fill(HIST("closureKstarRecoEP/ResTrueFT0C"), centrality, std::cos(2.0 * (recoPsiFT0C - trueCollision.eventPlaneAngle())));
        if (moment.mcClosure) {
          if (isEffSample) {
            histos.fill(HIST("moment/hResMCTrueRef"), centrality, std::cos(2.0 * (psi - psiTrue)));
          }
          if (isTestSample) {
            histos.fill(HIST("moment/hResMCTrueTest"), centrality, std::cos(2.0 * (psi - psiTrue)));
            histos.fill(HIST("moment/ResFT0CTPC"), centrality, std::cos(2.0 * (recoPsiFT0C - psiTPC)));
            histos.fill(HIST("moment/ResFT0CFT0A"), centrality, std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
            histos.fill(HIST("moment/ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
          }
        }
      }
      const int momSubsample = getMomentSubsample(recCollision.globalIndex());
      auto recTracksThisColl = recTracks.sliceBy(perCollision, recCollision.globalIndex());

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
            if (isEffSample) {
              histos.fill(HIST("closureKstarRecoEP/hSparseRecTrueEff"), kstar.M(), kstar.Pt(), cosRec, centrality, 1.0);
              if (moment.mcClosure && passCandWindow<false>(kstar)) {
                fillMoments<false>(kstar.M(), kstar.Pt(), centrality, kMomMCRef, momSubsample, cosRec, getRho00Weight(cosGen, moment.mcRefRho00.value));
              }
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
          if (moment.mcClosure && passCandWindow<false>(kstar)) {
            fillMoments<false>(kstar.M(), kstar.Pt(), centrality, kMomSEUnlike, momSubsample, cosRec, weight);
          }

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
        if (isEffSample) {
          histos.fill(HIST("closureKstarRecoEP/hSparseGenEff"), kstarGen.M(), kstarGen.Pt(), cosGen, centrality, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closureKstarRecoEP/hSparseGen"), kstarGen.M(), kstarGen.Pt(), cosGen, centrality, weight);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureRecoEP, "Process MC closure of cos(theta*) for K*, reconstructed Psi_FT0C: SE, true rec, gen", false);

  void processMCClosureMERecoEP(CollisionMCRecMixTableRecoEP const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    ROOT::Math::PxPyPzMVector kaon, pion, kstar;
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeMCClosureRecoEP binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}};
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
      const double psi = getSAEventAngle(recoPsiFT0C);
      {
        const double psiTPC = recCollision.psiTPC();
        const double psiFT0A = recCollision.psiFT0A();
        const double qFT0C = recCollision.qFT0C();
        const double qTPC = recCollision.qTPC();
        const double qFT0A = recCollision.qFT0A();
        histos.fill(HIST("closurePhiRecoEP/ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0CTPC"), centrality, std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0CFT0A"), centrality, std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closurePhiRecoEP/ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
        histos.fill(HIST("closurePhiRecoEP/ResTrueFT0C"), centrality, std::cos(2.0 * (recoPsiFT0C - trueCollision.eventPlaneAngle())));
        if (moment.mcClosure && !doprocessMCClosureRecoEP) {
          if (isEffSample) {
            histos.fill(HIST("moment/hResMCTrueRef"), centrality, std::cos(2.0 * (psi - psiTrue)));
          }
          if (isTestSample) {
            histos.fill(HIST("moment/hResMCTrueTest"), centrality, std::cos(2.0 * (psi - psiTrue)));
            histos.fill(HIST("moment/ResFT0CTPC"), centrality, std::cos(2.0 * (recoPsiFT0C - psiTPC)));
            histos.fill(HIST("moment/ResFT0CFT0A"), centrality, std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
            histos.fill(HIST("moment/ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
          }
        }
      }
      const int momSubsample = getMomentSubsample(recCollision.globalIndex());
      auto recTracksThisColl = recTracks.sliceBy(perCollision, recCollision.globalIndex());

      for (const auto& track1 : recTracksThisColl) {
        if (!(track1.signed1Pt() > cfgCutCharge.value)) {
          continue;
        }
        if (!selectionTrackPhi(track1) || !selectionPIDPhi(track1) || !selectionClosurePrimary(track1)) {
          continue;
        }
        for (const auto& track2 : recTracksThisColl) {
          if (!(track2.signed1Pt() < cfgCutCharge.value)) {
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
            if (isEffSample) {
              histos.fill(HIST("closurePhiRecoEP/hSparseRecTrueEff"), phi.M(), phi.Pt(), cosRec, centrality, 1.0);
              if (moment.mcClosure && passCandWindow<true>(phi)) {
                fillMoments<true>(phi.M(), phi.Pt(), centrality, kMomMCRef, momSubsample, cosRec, getRho00Weight(cosGen, moment.mcRefRho00.value));
              }
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
          if (moment.mcClosure && passCandWindow<true>(phi)) {
            fillMoments<true>(phi.M(), phi.Pt(), centrality, kMomSEUnlike, momSubsample, cosRec, weight);
          }

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
        if (isEffSample) {
          histos.fill(HIST("closurePhiRecoEP/hSparseGenEff"), phiGen.M(), phiGen.Pt(), cosGen, centrality, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closurePhiRecoEP/hSparseGen"), phiGen.M(), phiGen.Pt(), cosGen, centrality, weight);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureRecoEPPhi, "Process MC closure of cos(theta*) for phi(1020), reconstructed Psi_FT0C: SE, true rec, gen", false);

  void processMCClosureMERecoEPPhi(CollisionMCRecMixTableRecoEP const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    ROOT::Math::PxPyPzMVector kaonPlus, kaonMinus, phi;
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeMCClosureRecoEP binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}};
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

  struct MomDaughter {
    float px;
    float py;
    float pz;
    float tpcNSigma;
    float tofNSigma;
    bool hasTOF;
  };

  bool passDerivedPID(const MomDaughter& daughter)
  {
    if (std::abs(daughter.tpcNSigma) > candCfg.nSigmaTPC) {
      return false;
    }
    return !daughter.hasTOF || std::abs(daughter.tofNSigma) <= candCfg.nSigmaTOF;
  }

  template <bool IsPhi>
  bool passCandWindow(const ROOT::Math::PxPyPzMVector& mother)
  {
    const double mass = mother.M();
    if constexpr (IsPhi) {
      return mass >= candCfg.phiMassMin && mass <= candCfg.phiMassMax && std::abs(mother.Rapidity()) <= phiSA.confRapidity;
    } else {
      return mass >= candCfg.kstarMassMin && mass <= candCfg.kstarMassMax && std::abs(mother.Rapidity()) <= confRapidity;
    }
  }

  template <bool IsPhi>
  bool makeCandPair(float px1, float py1, float pz1, float px2, float py2, float pz2, ROOT::Math::PxPyPzMVector& mother)
  {
    ROOT::Math::PxPyPzMVector first(px1, py1, pz1, MassKa);
    ROOT::Math::PxPyPzMVector second(px2, py2, pz2, IsPhi ? MassKa : MassPi);
    mother = first + second;
    if constexpr (IsPhi) {
      return passCandWindow<true>(mother);
    } else {
      return passCandWindow<false>(mother) && !isPhiReflection(first, second);
    }
  }

  template <typename TCand>
  MomDaughter getPosDaughter(const TCand& cand)
  {
    return {cand.posPx(), cand.posPy(), cand.posPz(), cand.posTPCNSigma(), cand.posTOFNSigma(), cand.posHasTOF()};
  }

  template <typename TCand>
  MomDaughter getNegDaughter(const TCand& cand)
  {
    return {cand.negPx(), cand.negPy(), cand.negPz(), cand.negTPCNSigma(), cand.negTOFNSigma(), cand.negHasTOF()};
  }

  template <bool IsPhi, typename TCands, typename F>
  void forEachDerivedCand(const TCands& cands, const F& f)
  {
    ROOT::Math::PxPyPzMVector mother;
    for (const auto& cand : cands) {
      const auto pos = getPosDaughter(cand);
      const auto neg = getNegDaughter(cand);
      if (!passDerivedPID(pos) || !passDerivedPID(neg)) {
        continue;
      }
      if constexpr (IsPhi) {
        if (makeCandPair<true>(pos.px, pos.py, pos.pz, neg.px, neg.py, neg.pz, mother)) {
          f(mother, ROOT::Math::PxPyPzMVector(neg.px, neg.py, neg.pz, MassKa));
        }
      } else {
        const auto& kaon = cand.flavour() > 0 ? pos : neg;
        const auto& pion = cand.flavour() > 0 ? neg : pos;
        if (makeCandPair<false>(kaon.px, kaon.py, kaon.pz, pion.px, pion.py, pion.pz, mother)) {
          f(mother, ROOT::Math::PxPyPzMVector(kaon.px, kaon.py, kaon.pz, MassKa));
        }
      }
    }
  }

  void processMomentDerivedSE(aod::SpinAlignEvent const& event, aod::SAKstarCands const& kstars, aod::SAPhiCands const& phis)
  {
    const double psiSA = getSAEventAngle(event.psiFT0C());
    const double centrality = event.centFT0C();
    const int subsample = getMomentSubsample(event.globalIndex());
    histos.fill(HIST("moment/ResFT0CTPC"), centrality, std::cos(2.0 * (event.psiFT0C() - event.psiTPC())));
    histos.fill(HIST("moment/ResFT0CFT0A"), centrality, std::cos(2.0 * (event.psiFT0C() - event.psiFT0A())));
    histos.fill(HIST("moment/ResFT0ATPC"), centrality, std::cos(2.0 * (event.psiTPC() - event.psiFT0A())));
    forEachDerivedCand<false>(kstars, [&](const auto& mother, const auto& daughter) {
      fillMoments<false>(mother.M(), mother.Pt(), centrality, kMomSEUnlike, subsample, getCosThetaStar(mother, daughter, getSAAxis(mother, psiSA)));
    });
    forEachDerivedCand<true>(phis, [&](const auto& mother, const auto& daughter) {
      fillMoments<true>(mother.M(), mother.Pt(), centrality, kMomSEUnlike, subsample, getCosThetaStar(mother, daughter, getSAAxis(mother, psiSA)));
    });
  }
  PROCESS_SWITCH(Kstarpbpb, processMomentDerivedSE, "rho00 derived: SE moments (sample 0), K* and phi", false);

  struct ReplCand {
    float mass;
    float pt;
    float psi;
    float psiRP;
    float weight;
    int64_t bcId;
    double ux, uy, uz;
  };

  struct ReplCellBuffer {
    std::vector<ReplCand> cands;
    std::size_t next = 0;
  };

  struct ReplPool {
    std::unordered_map<int64_t, ReplCellBuffer> kstar, phi;
  };

  struct ReplEvent {
    int64_t index = -1;
    int64_t bcId = -1;
    double psi = 0.;
    double psiRP = 0.;
    bool hasRP = false;
    double centrality = 0.;
    std::vector<std::pair<int64_t, ReplCand>> kstar, phi;
  };

  using BinningTypeRepl = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  static constexpr int64_t ReplCellStride = 4096;
  static constexpr double ReplEtaOffset = 10.0;

  PresliceOptional<aod::SAKstarCands> kstarPerEventDerived = aod::savmcand::spinAlignEventId;
  PresliceOptional<aod::SAPhiCands> phiPerEventDerived = aod::savmcand::spinAlignEventId;

  int64_t getReplCell(double pt, double eta, double phi)
  {
    double phi02 = std::fmod(phi, o2::constants::math::TwoPI);
    if (phi02 < 0.) {
      phi02 += o2::constants::math::TwoPI;
    }
    const auto iPt = static_cast<int64_t>(std::floor(pt / moment.replPtWidth.value));
    const auto iEta = static_cast<int64_t>(std::floor((eta + ReplEtaOffset) / moment.replEtaWidth.value));
    const auto iPhi = static_cast<int64_t>(std::floor(phi02 / moment.replPhiWidth.value));
    return (iPt * ReplCellStride + iEta) * ReplCellStride + iPhi;
  }

  template <bool IsPhi>
  void addReplCand(ReplEvent& event, const ROOT::Math::PxPyPzMVector& mother, const ROOT::Math::PxPyPzMVector& daughter, double weight = 1.0)
  {
    if (!passCandWindow<IsPhi>(mother)) {
      return;
    }
    ROOT::Math::Boost boost{mother.BoostToCM()};
    const auto dir = boost(daughter).Vect().Unit();
    ReplCand cand{.mass = static_cast<float>(mother.M()), .pt = static_cast<float>(mother.Pt()), .psi = static_cast<float>(event.psi), .psiRP = static_cast<float>(event.psiRP), .weight = static_cast<float>(weight), .bcId = event.bcId, .ux = dir.X(), .uy = dir.Y(), .uz = dir.Z()};
    auto& cands = IsPhi ? event.phi : event.kstar;
    cands.emplace_back(getReplCell(mother.Pt(), mother.Eta(), mother.Phi()), cand);
  }

  template <bool IsPhi>
  void mixReplacement(const ReplEvent& event, ReplPool& pool)
  {
    const auto& targets = IsPhi ? event.phi : event.kstar;
    auto& cells = IsPhi ? pool.phi : pool.kstar;
    const auto targetAxis = getEPAxis(event.psi);
    const int subsample = getMomentSubsample(event.index);
    for (const auto& [cell, target] : targets) {
      const auto found = cells.find(cell);
      if (found == cells.end() || found->second.cands.empty()) {
        if constexpr (IsPhi) {
          histos.fill(HIST("moment/phi/hReplacement"), 1);
        } else {
          histos.fill(HIST("moment/kstar/hReplacement"), 1);
        }
        continue;
      }
      if constexpr (IsPhi) {
        histos.fill(HIST("moment/phi/hReplacement"), 0);
      } else {
        histos.fill(HIST("moment/kstar/hReplacement"), 0);
      }
      const auto& donors = found->second.cands;
      for (int iDraw = 0; iDraw < moment.replNDraws.value; ++iDraw) {
        const auto& donor = donors[rndGen.Integer(static_cast<unsigned int>(donors.size()))];
        if (event.bcId >= 0 && donor.bcId == event.bcId) {
          continue;
        }
        const auto axis = cfgSAFrame.value == kRandom3D ? getRandom3DAxis() : targetAxis;
        const double cosThetaStar = donor.ux * axis.X() + donor.uy * axis.Y() + donor.uz * axis.Z();
        const double weight = target.weight * donor.weight;
        const double cos2DeltaPsi = std::cos(2.0 * (event.psi - donor.psi));
        fillMoments<IsPhi>(donor.mass, donor.pt, event.centrality, kMomMERepl, subsample, cosThetaStar, weight);
        histos.fill(HIST("moment/hReplCos2DeltaPsi"), event.centrality, cos2DeltaPsi);
        if constexpr (IsPhi) {
          histos.fill(HIST("moment/phi/hXD"), donor.mass, donor.pt, event.centrality, kMomMERepl, subsample, weight * cos2DeltaPsi);
        } else {
          histos.fill(HIST("moment/kstar/hXD"), donor.mass, donor.pt, event.centrality, kMomMERepl, subsample, weight * cos2DeltaPsi);
        }
        if (event.hasRP) {
          const double cos2EPRP = std::cos(2.0 * (event.psi - donor.psiRP));
          if constexpr (IsPhi) {
            histos.fill(HIST("moment/phi/hXC"), donor.mass, donor.pt, event.centrality, kMomMERepl, subsample, weight * cos2EPRP);
          } else {
            histos.fill(HIST("moment/kstar/hXC"), donor.mass, donor.pt, event.centrality, kMomMERepl, subsample, weight * cos2EPRP);
          }
        }
      }
    }
    const auto capacity = static_cast<std::size_t>(moment.replCellCapacity.value);
    for (const auto& [cell, cand] : targets) {
      auto& buffer = cells[cell];
      if (buffer.cands.size() < capacity) {
        buffer.cands.push_back(cand);
      } else {
        buffer.cands[buffer.next] = cand;
        buffer.next = (buffer.next + 1) % capacity;
      }
    }
  }

  template <typename TTracks>
  void buildReplCandsKstarMC(ReplEvent& event, TTracks const& tracks)
  {
    ROOT::Math::PxPyPzMVector kaon, pion;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1) || !selectionClosurePrimary(track1)) {
        continue;
      }
      if ((!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) || (isTOFOnly && !selectionPID2(track1, 0))) {
        continue;
      }
      for (const auto& track2 : tracks) {
        if (track1.sign() * track2.sign() >= 0) {
          continue;
        }
        if (!selectionTrack(track2) || !selectionClosurePrimary(track2)) {
          continue;
        }
        if ((!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) || (isTOFOnly && !selectionPID2(track2, 1))) {
          continue;
        }
        kaon = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), MassKa);
        pion = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), MassPi);
        if (isPhiReflection(kaon, pion)) {
          continue;
        }
        double weight = 1.0;
        if (isTrueKstarPair(track1, track2)) {
          const auto mcKaon = track1.mcParticle();
          const auto mcPion = track2.mcParticle();
          const ROOT::Math::PxPyPzMVector kaonGen(mcKaon.px(), mcKaon.py(), mcKaon.pz(), MassKa);
          const ROOT::Math::PxPyPzMVector pionGen(mcPion.px(), mcPion.py(), mcPion.pz(), MassPi);
          const auto kstarGen = kaonGen + pionGen;
          weight = getInjectionWeight(getCosThetaStar(kstarGen, kaonGen, getSAAxis(kstarGen, event.psiRP)));
        }
        addReplCand<false>(event, kaon + pion, kaon, weight);
      }
    }
  }

  template <typename TTracks>
  void buildReplCandsPhiMC(ReplEvent& event, TTracks const& tracks)
  {
    ROOT::Math::PxPyPzMVector kaonPlus, kaonMinus;
    for (const auto& track1 : tracks) {
      if (!(track1.signed1Pt() > cfgCutCharge.value) || !selectionTrackPhi(track1) || !selectionPIDPhi(track1) || !selectionClosurePrimary(track1)) {
        continue;
      }
      for (const auto& track2 : tracks) {
        if (!(track2.signed1Pt() < cfgCutCharge.value) || !selectionTrackPhi(track2) || !selectionPIDPhi(track2) || !selectionClosurePrimary(track2)) {
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
        double weight = 1.0;
        if (isTruePhiPair(track1, track2)) {
          const auto mcPlus = track1.mcParticle();
          const auto mcMinus = track2.mcParticle();
          const ROOT::Math::PxPyPzMVector kaonPlusGen(mcPlus.px(), mcPlus.py(), mcPlus.pz(), MassKa);
          const ROOT::Math::PxPyPzMVector kaonMinusGen(mcMinus.px(), mcMinus.py(), mcMinus.pz(), MassKa);
          const auto phiGen = kaonPlusGen + kaonMinusGen;
          weight = getInjectionWeight(getCosThetaStar(phiGen, kaonMinusGen, getSAAxis(phiGen, event.psiRP)));
        }
        addReplCand<true>(event, kaonPlus + kaonMinus, kaonMinus, weight);
      }
    }
  }

  template <bool IsPhi>
  void runMCClosureReplacementME(CollisionMCRecMixTableRecoEP const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels)
  {
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());
    BinningTypeRepl binning{{moment.replAxisVertex, moment.axisCentrality}, true};
    std::vector<ReplPool> pools(binning.getAllBinsCount());
    for (const auto& collision : collisions) {
      if (!selectionEventMCRecoEP(collision) || !collision.has_mcCollision()) {
        continue;
      }
      const auto mcCollisionId = collision.mcCollisionId();
      if (!isClosureTestSample(mcCollisionId) || nRecPerMc[mcCollisionId] != 1) {
        continue;
      }
      const int bin = binning.getBin(std::make_tuple(collision.posZ(), collision.centFT0C()));
      if (bin < 0) {
        continue;
      }
      auto tracksThisEvent = tracks.sliceBy(perCollision, collision.globalIndex());
      ReplEvent event;
      event.index = collision.globalIndex();
      event.bcId = collision.bcId();
      event.psi = getSAEventAngle(collision.psiFT0C());
      event.psiRP = getSAEventAngle(collision.mcCollision_as<aod::McCollisions>().eventPlaneAngle());
      event.hasRP = true;
      event.centrality = collision.centFT0C();
      if constexpr (IsPhi) {
        buildReplCandsPhiMC(event, tracksThisEvent);
      } else {
        buildReplCandsKstarMC(event, tracksThisEvent);
      }
      mixReplacement<IsPhi>(event, pools[bin]);
    }
  }

  void processMCClosureReplacementME(CollisionMCRecMixTableRecoEP const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    runMCClosureReplacementME<false>(collisions, tracks, mcCollisions, mcLabels);
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureReplacementME, "rho00 K* MC closure: replacement ME of the test sample with injected weights (sample 1), plane correlation hXD and hXC", false);

  void processMCClosureReplacementMEPhi(CollisionMCRecMixTableRecoEP const& collisions, FilTrackMCRecTable const& tracks, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    runMCClosureReplacementME<true>(collisions, tracks, mcCollisions, mcLabels);
  }
  PROCESS_SWITCH(Kstarpbpb, processMCClosureReplacementMEPhi, "rho00 phi MC closure: replacement ME of the test sample with injected weights (sample 1), plane correlation hXD and hXC", false);

  void processMomentDerivedMERepl(aod::SpinAlignEvents const& events, aod::SAKstarCands const& kstars, aod::SAPhiCands const& phis)
  {
    BinningTypeRepl binning{{moment.replAxisVertex, moment.axisCentrality}, true};
    std::vector<ReplPool> pools(binning.getAllBinsCount());
    for (const auto& ev : events) {
      const int bin = binning.getBin(std::make_tuple(ev.posZ(), ev.centFT0C()));
      if (bin < 0) {
        continue;
      }
      ReplEvent event;
      event.index = ev.globalIndex();
      event.psi = getSAEventAngle(ev.psiFT0C());
      event.centrality = ev.centFT0C();
      forEachDerivedCand<false>(kstars.sliceBy(kstarPerEventDerived, ev.globalIndex()), [&](const auto& mother, const auto& daughter) {
        addReplCand<false>(event, mother, daughter);
      });
      forEachDerivedCand<true>(phis.sliceBy(phiPerEventDerived, ev.globalIndex()), [&](const auto& mother, const auto& daughter) {
        addReplCand<true>(event, mother, daughter);
      });
      mixReplacement<false>(event, pools[bin]);
      mixReplacement<true>(event, pools[bin]);
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMomentDerivedMERepl, "rho00 derived: ME reference by candidate replacement (sample 1), K* and phi", false);

  struct V2Track {
    float px;
    float py;
    float pz;
    float phi;
    int8_t sign;
    bool isKaon;
    bool isPion;
  };

  struct V2Event {
    int64_t index = -1;
    double psi = 0.;
    double q = 0.;
    double centrality = 0.;
    std::vector<V2Track> tracks;
  };

  using BinningTypeV2Derived = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  PresliceOptional<aod::SAVMTracks> trackPerEventV2 = aod::savmtrack::spinAlignEventId;

  template <typename TEvent, typename TTracks>
  V2Event buildV2Event(const TEvent& event, const TTracks& tracks)
  {
    V2Event v2Event;
    v2Event.index = event.globalIndex();
    v2Event.psi = event.psiFT0C();
    v2Event.q = event.qFT0C();
    v2Event.centrality = event.centFT0C();
    for (const auto& track : tracks) {
      const bool kaon = track.isKaon() && passDerivedPID(MomDaughter{track.px(), track.py(), track.pz(), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.hasTOF()});
      const bool pion = track.isPion() && passDerivedPID(MomDaughter{track.px(), track.py(), track.pz(), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.hasTOF()});
      if (!kaon && !pion) {
        continue;
      }
      v2Event.tracks.push_back({track.px(), track.py(), track.pz(), static_cast<float>(std::atan2(track.py(), track.px())), track.sign(), kaon, pion});
    }
    return v2Event;
  }

  void fillV2LikeSignSE(const V2Event& event)
  {
    ROOT::Math::PxPyPzMVector mother;
    const auto& tracks = event.tracks;
    for (std::size_t i = 0; i < tracks.size(); ++i) {
      for (std::size_t j = 0; j < tracks.size(); ++j) {
        if (i == j || tracks[i].sign != tracks[j].sign) {
          continue;
        }
        if (tracks[i].isKaon && tracks[j].isPion && makeCandPair<false>(tracks[i].px, tracks[i].py, tracks[i].pz, tracks[j].px, tracks[j].py, tracks[j].pz, mother)) {
          if (tracks[i].sign > 0) {
            histos.fill(HIST("v2derived/kstar/hSELikePP"), mother.M(), mother.Pt(), event.centrality);
          } else {
            histos.fill(HIST("v2derived/kstar/hSELikeNN"), mother.M(), mother.Pt(), event.centrality);
          }
        }
        if (j > i && tracks[i].isKaon && tracks[j].isKaon && makeCandPair<true>(tracks[i].px, tracks[i].py, tracks[i].pz, tracks[j].px, tracks[j].py, tracks[j].pz, mother)) {
          if (tracks[i].sign > 0) {
            histos.fill(HIST("v2derived/phi/hSELikePP"), mother.M(), mother.Pt(), event.centrality);
          } else {
            histos.fill(HIST("v2derived/phi/hSELikeNN"), mother.M(), mother.Pt(), event.centrality);
          }
        }
      }
    }
  }

  static double v2DaughterTerm(const V2Event& event, const V2Track& track, double phiMother)
  {
    return event.q * std::cos(2.0 * (track.phi - event.psi)) * std::cos(2.0 * (track.phi - phiMother));
  }

  void fillV2MixedKstar(const V2Event& a, const V2Event& b, double centrality)
  {
    ROOT::Math::PxPyPzMVector mother;
    for (const auto& kaon : a.tracks) {
      if (!kaon.isKaon) {
        continue;
      }
      for (const auto& pion : b.tracks) {
        if (!pion.isPion || !makeCandPair<false>(kaon.px, kaon.py, kaon.pz, pion.px, pion.py, pion.pz, mother)) {
          continue;
        }
        if (kaon.sign != pion.sign) {
          const double term = v2DaughterTerm(a, kaon, mother.Phi()) + v2DaughterTerm(b, pion, mother.Phi());
          histos.fill(HIST("v2derived/kstar/hMEUnlike"), mother.M(), mother.Pt(), term, centrality);
        } else if (kaon.sign > 0) {
          histos.fill(HIST("v2derived/kstar/hMELikePP"), mother.M(), mother.Pt(), centrality);
        } else {
          histos.fill(HIST("v2derived/kstar/hMELikeNN"), mother.M(), mother.Pt(), centrality);
        }
      }
    }
  }

  void fillV2MixedPhi(const V2Event& a, const V2Event& b, double centrality)
  {
    ROOT::Math::PxPyPzMVector mother;
    for (const auto& kaonA : a.tracks) {
      if (!kaonA.isKaon) {
        continue;
      }
      for (const auto& kaonB : b.tracks) {
        if (!kaonB.isKaon || !makeCandPair<true>(kaonA.px, kaonA.py, kaonA.pz, kaonB.px, kaonB.py, kaonB.pz, mother)) {
          continue;
        }
        if (kaonA.sign != kaonB.sign) {
          const double term = v2DaughterTerm(a, kaonA, mother.Phi()) + v2DaughterTerm(b, kaonB, mother.Phi());
          histos.fill(HIST("v2derived/phi/hMEUnlike"), mother.M(), mother.Pt(), term, centrality);
        } else if (kaonA.sign > 0) {
          histos.fill(HIST("v2derived/phi/hMELikePP"), mother.M(), mother.Pt(), centrality);
        } else {
          histos.fill(HIST("v2derived/phi/hMELikeNN"), mother.M(), mother.Pt(), centrality);
        }
      }
    }
  }

  void processV2DerivedSE(aod::SpinAlignEvent const& event, aod::SAKstarCands const& kstars, aod::SAPhiCands const& phis, aod::SAVMTracks const& tracks)
  {
    const double centrality = event.centFT0C();
    const double psi = event.psiFT0C();
    const double q = event.qFT0C();
    histos.fill(HIST("v2derived/ResFT0CTPCSP"), centrality, event.qFT0C() * event.qTPC() * std::cos(2.0 * (event.psiFT0C() - event.psiTPC())));
    histos.fill(HIST("v2derived/ResFT0CFT0ASP"), centrality, event.qFT0C() * event.qFT0A() * std::cos(2.0 * (event.psiFT0C() - event.psiFT0A())));
    histos.fill(HIST("v2derived/ResFT0ATPCSP"), centrality, event.qFT0A() * event.qTPC() * std::cos(2.0 * (event.psiFT0A() - event.psiTPC())));

    forEachDerivedCand<false>(kstars, [&](const auto& mother, const auto&) {
      histos.fill(HIST("v2derived/kstar/hSEUnlike"), mother.M(), mother.Pt(), q * std::cos(2.0 * (mother.Phi() - psi)), centrality);
    });
    forEachDerivedCand<true>(phis, [&](const auto& mother, const auto&) {
      histos.fill(HIST("v2derived/phi/hSEUnlike"), mother.M(), mother.Pt(), q * std::cos(2.0 * (mother.Phi() - psi)), centrality);
    });

    fillV2LikeSignSE(buildV2Event(event, tracks));
  }
  PROCESS_SWITCH(Kstarpbpb, processV2DerivedSE, "v2 derived: same event (u2 Q2* vs mass, like-sign pairs, SP resolution), K* and phi", false);

  void processV2DerivedME(aod::SpinAlignEvents const& events, aod::SAVMTracks const& tracks)
  {
    BinningTypeV2Derived binning{{v2Derived.mixAxisVertex, v2Derived.axisCentrality}, true};
    std::vector<std::deque<V2Event>> pools(binning.getAllBinsCount());
    for (const auto& event : events) {
      const int bin = binning.getBin(std::make_tuple(event.posZ(), event.centFT0C()));
      if (bin < 0) {
        continue;
      }
      auto tracksThisEvent = tracks.sliceBy(trackPerEventV2, event.globalIndex());
      V2Event current = buildV2Event(event, tracksThisEvent);
      for (const auto& previous : pools[bin]) {
        histos.fill(HIST("v2derived/hMixedEventPairs"), current.centrality);
        fillV2MixedKstar(current, previous, current.centrality);
        fillV2MixedKstar(previous, current, current.centrality);
        fillV2MixedPhi(current, previous, current.centrality);
      }
      pools[bin].push_back(std::move(current));
      if (static_cast<int>(pools[bin].size()) > v2Derived.nMixedEvents.value) {
        pools[bin].pop_front();
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processV2DerivedME, "v2 derived: mixed event (background v2, Eq. 5; like-sign normalisation, Eq. 6), K* and phi", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Kstarpbpb>(cfgc)};
}
