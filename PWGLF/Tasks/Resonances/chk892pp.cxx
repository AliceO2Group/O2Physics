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

/// \file chk892pp.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <Bong-Hwi.Lim@cern.ch>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>
// #include <TDatabasePDG.h> // FIXME
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
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
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/RotationZ.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TParticlePDG.h"
#include "TRandom3.h"
#include "TVector2.h"
#include <TMath.h>
#include <TPDGCode.h> // FIXME

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::aod::rctsel;

namespace
{
template <typename V0T>
inline bool getTruthK0sAndGenKinematics(V0T const& v0, double& ptgen, double& ygen)
{
  if (!v0.has_mcParticle())
    return false;
  auto mcPart = v0.template mcParticle_as<aod::McParticles>();
  if (mcPart.pdgCode() != kK0Short)
    return false;
  ptgen = mcPart.pt();
  ygen = mcPart.y();
  return true;
}
} // namespace

struct Chk892pp {
  enum BinType : unsigned int {
    kKstarP = 0,
    kKstarN,
    kKstarP_Mix,
    kKstarN_Mix,
    kKstarP_Rot,
    kKstarN_Rot,
    kTYEnd
  };

  enum EvtStep {
    kAll = 0,
    kZvtx,
    kINELgt0,
    kAssocReco,
    kNSteps
  };

  const int nSteps = static_cast<int>(EvtStep::kNSteps);

  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<aod::V0s> perCollisionV0 = aod::v0data::collisionId;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::PVMults>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using V0Candidates = aod::V0Datas;

  // for MC reco
  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Join<TrackCandidates, aod::McTrackLabels>; //, aod::McParticles>;
  using MCV0Candidates = soa::Join<V0Candidates, aod::McV0Labels>;

  // for MC truth
  using MCTrueEventCandidates = aod::McCollisions;
  using MCTrueTrackCandidates = aod::McParticles;

  using LorentzVectorSetXYZM = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::ccdb::CcdbApi ccdbApi;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  } CCDBConfig;
  // Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  // Configurables
  struct : ConfigurableGroup {
    ConfigurableAxis cfgBinsPt{"cfgBinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsPtQA{"cfgBinsPtQA", {VARIABLE_WIDTH, 0.0, 0.8, 1.3, 1.8, 2.3, 2.8, 3.4, 4.0, 5.0, 6.0, 7.0, 8.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsCent{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
    ConfigurableAxis cfgBinsVtxZ{"cfgBinsVtxZ", {VARIABLE_WIDTH, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "Binning of the z-vertex axis"};
    Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};
    Configurable<int> cNbinsDivQA{"cNbinsDivQA", 1, "Integer to divide the number of bins for QA"};
  } AxisConfig;

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  struct : ConfigurableGroup {
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<int> cfgEvtOccupancyInTimeRangeMax{"cfgEvtOccupancyInTimeRangeMax", -1, "Evt sel: maximum track occupancy"};
    Configurable<int> cfgEvtOccupancyInTimeRangeMin{"cfgEvtOccupancyInTimeRangeMin", -1, "Evt sel: minimum track occupancy"};
    Configurable<bool> cfgEvtTriggerCheck{"cfgEvtTriggerCheck", false, "Evt sel: check for trigger"};
    Configurable<bool> cfgEvtOfflineCheck{"cfgEvtOfflineCheck", false, "Evt sel: check for offline selection"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", true, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtTFBorderCut{"cfgEvtTFBorderCut", true, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtUseITSTPCvertex{"cfgEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
    Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", false, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", false, "Evt sel: apply pileup rejection"};
    Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", true, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgincludeCentralityMC{"cfgincludeCentralityMC", false, "Include centrality in MC"};
    Configurable<bool> cfgEvtCollInTimeRangeStandard{"cfgEvtCollInTimeRangeStandard", false, "Evt sel: apply NoCollInTimeRangeStandard"};
    Configurable<float> cfgEventCentralityMin{"cfgEventCentralityMin", 0.0f, "Event sel: minimum centrality"};
    Configurable<float> cfgEventCentralityMax{"cfgEventCentralityMax", 100.0f, "Event sel: maximum centrality"};
    Configurable<bool> cfgEvtUseRCTFlagChecker{"cfgEvtUseRCTFlagChecker", false, "Evt sel: use RCT flag checker"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } EventCuts;
  RCTFlagsChecker rctChecker;

  Configurable<bool> cfgFillQAPlots{"cfgFillQAPlots", false, "Fill QA plots"};
  Configurable<int> cfgCentEst{"cfgCentEst", 2, "Centrality estimator, 1: FT0C, 2: FT0M"};

  /// PID Selections, pion
  struct : ConfigurableGroup {
    Configurable<bool> cfgTPConly{"cfgTPConly", false, "Use only TPC for PID"};                                     // bool
    Configurable<float> cfgMaxTPCnSigmaPion{"cfgMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};                 // TPC
    Configurable<float> cfgMaxTOFnSigmaPion{"cfgMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                 // TOF
    Configurable<float> cfgNsigmaCutCombinedPion{"cfgNsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"}; // Combined
    Configurable<bool> cfgTOFVeto{"cfgTOFVeto", false, "TOF Veto, if false, TOF is nessessary for PID selection"};  // TOF Veto
    Configurable<float> cfgTOFMinPt{"cfgTOFMinPt", 0.6, "Minimum TOF pT cut for Pion"};                             // TOF pT cut
  } PIDCuts;

  // Track selections
  struct : ConfigurableGroup {
    Configurable<float> cfgMinPtcut{"cfgMinPtcut", 0.15, "Track minium pt cut"};
    Configurable<float> cfgMaxEtacut{"cfgMaxEtacut", 0.8, "Track maximum eta cut"};
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor

    Configurable<bool> cfgpTdepDCAxyCut{"cfgpTdepDCAxyCut", true, "pT-dependent DCAxy cut"};
    Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
    Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};
    Configurable<float> cfgRatioTPCRowsOverFindableCls{"cfgRatioTPCRowsOverFindableCls", 0.0f, "TPC Crossed Rows to Findable Clusters"};
    Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 999.0, "ITS Chi2/NCl"};
    Configurable<float> cfgTPCChi2NCl{"cfgTPCChi2NCl", 999.0, "TPC Chi2/NCl"};
    Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
    Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
    Configurable<bool> cfgHasITS{"cfgHasITS", false, "Require ITS"};
    Configurable<bool> cfgHasTPC{"cfgHasTPC", false, "Require TPC"};
    Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
    // DCAr to PV
    Configurable<float> cfgMaxbDCArToPVcut{"cfgMaxbDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
    // DCAz to PV
    Configurable<float> cfgMaxbDCAzToPVcut{"cfgMaxbDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
  } TrackCuts;

  // Secondary Selection
  struct : ConfigurableGroup {
    Configurable<bool> cfgReturnFlag{"cfgReturnFlag", false, "Return Flag for debugging"};
    Configurable<bool> cfgSecondaryRequire{"cfgSecondaryRequire", true, "Secondary cuts on/off"};
    Configurable<bool> cfgSecondaryArmenterosCut{"cfgSecondaryArmenterosCut", true, "cut on Armenteros-Podolanski graph"};
    Configurable<bool> cfgSecondaryCrossMassHypothesisCut{"cfgSecondaryCrossMassHypothesisCut", false, "Apply cut based on the lambda mass hypothesis"};

    Configurable<bool> cfgByPassDauPIDSelection{"cfgByPassDauPIDSelection", true, "Bypass Daughters PID selection"};
    Configurable<float> cfgSecondaryDauDCAMax{"cfgSecondaryDauDCAMax", 1., "Maximum DCA Secondary daughters to PV"};
    Configurable<float> cfgSecondaryDauPosDCAtoPVMin{"cfgSecondaryDauPosDCAtoPVMin", 0.1, "Minimum DCA Secondary positive daughters to PV"};
    Configurable<float> cfgSecondaryDauNegDCAtoPVMin{"cfgSecondaryDauNegDCAtoPVMin", 0.1, "Minimum DCA Secondary negative daughters to PV"};

    Configurable<float> cfgSecondaryPtMin{"cfgSecondaryPtMin", 0.f, "Minimum transverse momentum of Secondary"};
    Configurable<float> cfgSecondaryRapidityMax{"cfgSecondaryRapidityMax", 0.8, "Maximum rapidity of Secondary"};
    Configurable<float> cfgSecondaryRadiusMin{"cfgSecondaryRadiusMin", 0, "Minimum transverse radius of Secondary"};
    Configurable<float> cfgSecondaryRadiusMax{"cfgSecondaryRadiusMax", 999.9, "Maximum transverse radius of Secondary"};
    Configurable<float> cfgSecondaryCosPAMin{"cfgSecondaryCosPAMin", 0.995, "Mininum cosine pointing angle of Secondary"};
    Configurable<float> cfgSecondaryDCAtoPVMax{"cfgSecondaryDCAtoPVMax", 0.4, "Maximum DCA Secondary to PV"};
    Configurable<float> cfgSecondaryProperLifetimeMax{"cfgSecondaryProperLifetimeMax", 20, "Maximum Secondary Lifetime"};
    Configurable<float> cfgSecondaryparamArmenterosCut{"cfgSecondaryparamArmenterosCut", 0.2, "parameter for Armenteros Cut"};
    Configurable<float> cfgSecondaryMassWindow{"cfgSecondaryMassWindow", 0.03, "Secondary inv mass selciton window"};
    Configurable<float> cfgSecondaryCrossMassCutWindow{"cfgSecondaryCrossMassCutWindow", 0.05, "Secondary inv mass selection window with (anti)lambda hypothesis"};
  } SecondaryCuts;

  // K* selection
  struct : ConfigurableGroup {
    Configurable<float> cfgKstarMaxRap{"cfgKstarMaxRap", 0.5, "Kstar maximum rapidity"};
    Configurable<float> cfgKstarMinRap{"cfgKstarMinRap", -0.5, "Kstar minimum rapidity"};
  } KstarCuts;

  // Bkg estimation
  struct : ConfigurableGroup {
    Configurable<bool> cfgFillRotBkg{"cfgFillRotBkg", true, "Fill rotated background"};
    Configurable<float> cfgMinRot{"cfgMinRot", 5.0 * constants::math::PI / 6.0, "Minimum of rotation"};
    Configurable<float> cfgMaxRot{"cfgMaxRot", 7.0 * constants::math::PI / 6.0, "Maximum of rotation"};
    Configurable<bool> cfgRotPion{"cfgRotPion", true, "Rotate pion"};
    Configurable<int> cfgNrotBkg{"cfgNrotBkg", 4, "Number of rotated copies (background) per each original candidate"};
  } BkgEstimationConfig;

  Configurable<bool> cfgTruthUseInelGt0{"cfgTruthUseInelGt0", true, "Truth denominator: require INEL>0"};
  Configurable<bool> cfgTruthIncludeZvtx{"cfgTruthIncludeZvtx", true, "Truth denominator: also require |vtxz|<cfgEvtZvtx"};

  float lCentrality;

  // PDG code
  int kPDGK0s = kK0Short;
  int kPDGK0 = kK0;
  int kKstarPlus = o2::constants::physics::Pdg::kKPlusStar892;
  // int kPiPlus = 211;

  void init(o2::framework::InitContext&)
  {
    lCentrality = -999;

    // setCuts(float zvtxMax, bool checkTrigger, bool checkOffline, bool checkRun3, bool triggerTVXsel = false, int trackOccupancyInTimeRangeMax = -1, int trackOccupancyInTimeRangeMin = -1)
    colCuts.setCuts(EventCuts.cfgEvtZvtx, EventCuts.cfgEvtTriggerCheck, EventCuts.cfgEvtOfflineCheck, /*checkRun3*/ true, EventCuts.cfgEvtTriggerTVXSel, EventCuts.cfgEvtOccupancyInTimeRangeMax, EventCuts.cfgEvtOccupancyInTimeRangeMin);
    colCuts.init(&histos);
    colCuts.setTriggerTVX(EventCuts.cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(EventCuts.cfgEvtTFBorderCut);
    colCuts.setApplyNoITSROBorderCut(EventCuts.cfgEvtNoITSROBorderCut);
    colCuts.setApplyITSTPCvertex(EventCuts.cfgEvtUseITSTPCvertex);
    colCuts.setApplyZvertexTimedifference(EventCuts.cfgEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(EventCuts.cfgEvtPileupRejection);
    colCuts.setApplyCollInTimeRangeStandard(EventCuts.cfgEvtCollInTimeRangeStandard);
    colCuts.printCuts();

    rctChecker.init(EventCuts.cfgEvtRCTFlagCheckerLabel, EventCuts.cfgEvtRCTFlagCheckerZDCCheck, EventCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    AxisSpec centAxis = {AxisConfig.cfgBinsCent, "T0M (%)"};
    AxisSpec vtxzAxis = {AxisConfig.cfgBinsVtxZ, "Z Vertex (cm)"};
    AxisSpec epAxis = {100, -1.0 * constants::math::PI, constants::math::PI};
    AxisSpec epresAxis = {100, -1.02, 1.02};
    AxisSpec ptAxis = {AxisConfig.cfgBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {AxisConfig.cfgBinsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec radiusAxis = {50, 0, 5, "Radius (cm)"};
    AxisSpec cpaAxis = {50, 0.95, 1.0, "CPA"};
    AxisSpec tauAxis = {250, 0, 25, "Lifetime (cm)"};
    AxisSpec dcaAxis = {200, 0, 2, "DCA (cm)"};
    AxisSpec dcaxyAxis = {200, 0, 2, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {200, 0, 2, "DCA_{#it{z}} (cm)"};
    AxisSpec yAxis = {100, -1, 1, "Rapidity"};
    AxisSpec invMassAxisK0s = {800 / AxisConfig.cNbinsDiv, 0.46, 0.54, "Invariant Mass (GeV/#it{c}^2)"};  // K0s ~497.611
    AxisSpec invMassAxisReso = {900 / AxisConfig.cNbinsDiv, 0.5f, 1.4f, "Invariant Mass (GeV/#it{c}^2)"}; // chK(892) ~892
    AxisSpec pidQAAxis = {130 / AxisConfig.cNbinsDivQA, -6.5, 6.5};
    AxisSpec dataTypeAxis = {9, 0, 9, "Histogram types"};
    AxisSpec mcTypeAxis = {4, 0, 4, "Histogram types"};
    AxisSpec cutAxis{14, 0.5, 14.5, "Check"}; // 1..14

    // THnSparse
    AxisSpec axisType = {BinType::kTYEnd, 0, BinType::kTYEnd, "Type of bin with charge and mix"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};

    histos.add("QA/K0sCutCheck", "Check K0s cut", HistType::kTH1D, {AxisSpec{13, -0.5, 12.5, "Check"}});
    histos.add("QA/K0sCutFlow", "Check K0s cut (first-fail or pass)", HistType::kTH1F, {cutAxis});
    auto hcut = histos.get<TH1>(HIST("QA/K0sCutFlow"));
    hcut->GetXaxis()->SetBinLabel(1, "TOTAL");
    hcut->GetXaxis()->SetBinLabel(2, "PASS");
    hcut->GetXaxis()->SetBinLabel(3, "DauDCA>max");
    hcut->GetXaxis()->SetBinLabel(4, "PosDCAtoPV<min");
    hcut->GetXaxis()->SetBinLabel(5, "NegDCAtoPV<min");
    hcut->GetXaxis()->SetBinLabel(6, "pT<min");
    hcut->GetXaxis()->SetBinLabel(7, "|y|>max");
    hcut->GetXaxis()->SetBinLabel(8, "R<min or R>max");
    hcut->GetXaxis()->SetBinLabel(9, "DCAtoPV>max");
    hcut->GetXaxis()->SetBinLabel(10, "cosPA<min");
    hcut->GetXaxis()->SetBinLabel(11, "ctau>max");
    hcut->GetXaxis()->SetBinLabel(12, "qtarm<a|alpha|");
    hcut->GetXaxis()->SetBinLabel(13, "|M(K0s)-m0|>win");
    hcut->GetXaxis()->SetBinLabel(14, "cross-mass veto");

    histos.add("QA/before/CentDist", "Centrality distribution", {HistType::kTH1D, {centAxis}});
    histos.add("QA/before/VtxZ", "Centrality distribution", {HistType::kTH1D, {vtxzAxis}});
    histos.add("QA/before/hEvent", "Number of Events", HistType::kTH1F, {{1, 0.5, 1.5}});

    if (BkgEstimationConfig.cfgFillRotBkg) {
      histos.add("QA/RotBkg/hRotBkg", "Rotated angle of rotated background", HistType::kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});
    }

    // Bachelor pion
    histos.add("QA/before/trkbpionDCAxy", "DCAxy distribution of bachelor pion candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/before/trkbpionDCAz", "DCAz distribution of bachelor pion candidates", HistType::kTH1D, {dcazAxis});
    histos.add("QA/before/trkbpionpT", "pT distribution of bachelor pion candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/before/trkbpionTPCPID", "TPC PID of bachelor pion candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trkbpionTOFPID", "TOF PID of bachelor pion candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trkbpionTPCTOFPID", "TPC-TOF PID map of bachelor pion candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});

    histos.add("QA/after/trkbpionDCAxy", "DCAxy distribution of bachelor pion candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/after/trkbpionDCAz", "DCAz distribution of bachelor pion candidates", HistType::kTH1D, {dcazAxis});
    histos.add("QA/after/trkbpionpT", "pT distribution of bachelor pion candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/after/trkbpionTPCPID", "TPC PID of bachelor pion candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trkbpionTOFPID", "TOF PID of bachelor pion candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trkbpionTPCTOFPID", "TPC-TOF PID map of bachelor pion candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});

    // Secondary pion 1
    histos.add("QA/before/trkppionTPCPID", "TPC PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trkppionTOFPID", "TOF PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trkppionTPCTOFPID", "TPC-TOF PID map of secondary pion 1 (positive) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});
    histos.add("QA/before/trkppionpT", "pT distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/before/trkppionDCAxy", "DCAxy distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/before/trkppionDCAz", "DCAz distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcazAxis});

    histos.add("QA/after/trkppionTPCPID", "TPC PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trkppionTOFPID", "TOF PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trkppionTPCTOFPID", "TPC-TOF PID map of secondary pion 1 (positive) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});
    histos.add("QA/after/trkppionpT", "pT distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/after/trkppionDCAxy", "DCAxy distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/after/trkppionDCAz", "DCAz distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcazAxis});

    // Secondary pion 2
    histos.add("QA/before/trknpionTPCPID", "TPC PID of secondary pion 2 (negative) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trknpionTOFPID", "TOF PID of secondary pion 2 (negative) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trknpionTPCTOFPID", "TPC-TOF PID map of secondary pion 2 (negative) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});
    histos.add("QA/before/trknpionpT", "pT distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/before/trknpionDCAxy", "DCAxy distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/before/trknpionDCAz", "DCAz distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {dcazAxis});

    histos.add("QA/after/trknpionTPCPID", "TPC PID of secondary pion 2 (negative) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trknpionTOFPID", "TOF PID of secondary pion 2 (negative) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trknpionTPCTOFPID", "TPC-TOF PID map of secondary pion 2 (negative) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});
    histos.add("QA/after/trknpionpT", "pT distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/after/trknpionDCAxy", "DCAxy distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/after/trknpionDCAz", "DCAz distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {dcazAxis});

    // K0s
    histos.add("QA/before/hDauDCASecondary", "DCA of daughters of secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/before/hDauPosDCAtoPVSecondary", "Pos DCA to PV of daughters secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/before/hDauNegDCAtoPVSecondary", "Neg DCA to PV of daughters secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/before/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/before/hy_Secondary", "Rapidity distribution of secondary resonance", HistType::kTH1D, {yAxis});
    histos.add("QA/before/hRadiusSecondary", "Radius distribution of secondary resonance", HistType::kTH1D, {radiusAxis});
    histos.add("QA/before/hCPASecondary", "Cosine pointing angle distribution of secondary resonance", HistType::kTH1D, {cpaAxis});
    histos.add("QA/before/hDCAtoPVSecondary", "DCA to PV distribution of secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/before/hPropTauSecondary", "Proper Lifetime distribution of secondary resonance", HistType::kTH1D, {tauAxis});
    histos.add("QA/before/hPtAsymSecondary", "pT asymmetry distribution of secondary resonance", HistType::kTH1D, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QA/before/hInvmassSecondary", "Invariant mass of unlike-sign secondary resonance", HistType::kTH1D, {invMassAxisK0s});

    histos.add("QA/after/hDauDCASecondary", "DCA of daughters of secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/after/hDauPosDCAtoPVSecondary", "Pos DCA to PV of daughters secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/after/hDauNegDCAtoPVSecondary", "Neg DCA to PV of daughters secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/after/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/after/hy_Secondary", "Rapidity distribution of secondary resonance", HistType::kTH1D, {yAxis});
    histos.add("QA/after/hRadiusSecondary", "Radius distribution of secondary resonance", HistType::kTH1D, {radiusAxis});
    histos.add("QA/after/hCPASecondary", "Cosine pointing angle distribution of secondary resonance", HistType::kTH1D, {cpaAxis});
    histos.add("QA/after/hDCAtoPVSecondary", "DCA to PV distribution of secondary resonance", HistType::kTH1D, {dcaAxis});
    histos.add("QA/after/hPropTauSecondary", "Proper Lifetime distribution of secondary resonance", HistType::kTH1D, {tauAxis});
    histos.add("QA/after/hPtAsymSecondary", "pT asymmetry distribution of secondary resonance", HistType::kTH1D, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QA/after/hInvmassSecondary", "Invariant mass of unlike-sign secondary resonance", HistType::kTH1D, {invMassAxisK0s});

    // Kstar
    // Invariant mass nSparse
    histos.add("QA/before/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
    histos.add("hInvmass_Kstar", "Invariant mass of unlike-sign chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_Kstar_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_K0s", "Invariant mass of unlike-sign K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s});

    // Mass QA (quick check)
    histos.add("QA/before/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
    histos.add("QA/before/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

    histos.add("QA/after/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
    histos.add("QA/after/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
    histos.add("QA/after/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

    // MC
    if (doprocessMC) {

      histos.add("QACent_woCut", "Centrality without cut", HistType::kTH1F, {centAxis});
      histos.add("QACent_woCentCut", "Centrality without cent cut", HistType::kTH1F, {centAxis});
      histos.add("QACent_wCentCut", "Centrality with cent cut", HistType::kTH1F, {centAxis});
      histos.add("QAvtxz_woCut", "z-vertex without cut", HistType::kTH1F, {vtxzAxis});
      histos.add("QAvtxz_wVtxzCut", "z-vertex with vtxz cut", HistType::kTH1F, {vtxzAxis});

      histos.add("EffK0s/genK0s", "Gen K0s (|y<0.8|)", HistType::kTH2F, {ptAxisQA, centAxis});
      histos.add("EffK0s/recoK0s", "Reco K0s (|y<0.8|)", HistType::kTH2F, {ptAxisQA, centAxis});

      histos.add("EffKstar/genKstar", "Gen Kstar (|y|<0.5)", HistType::kTH2F, {ptAxisQA, centAxis});
      histos.add("EffKstar/recoKstar", "Kstar Reco matched (final all)", HistType::kTH2F, {ptAxisQA, centAxis});

      histos.add("Correction/sigLoss_den", "Gen Kstar (|y|<0.5) in truth class", HistType::kTH2F, {ptAxisQA, centAxis});
      histos.add("Correction/sigLoss_num", "Gen Kstar (|y|<0.5, selected events) in reco class", HistType::kTH2F, {ptAxisQA, centAxis});
      histos.add("Correction/EF_den", "Gen events (truth class)", HistType::kTH1F, {centAxis});
      histos.add("Correction/EF_num", "Reco events (selected events)", HistType::kTH1F, {centAxis});
      histos.add("Correction/MCTruthCent_all", "MC truth FT0M centrality (all mcCollisions)", HistType::kTH1F, {centAxis});
      histos.add("Correction/MCTruthCent_cut", "MC truth FT0M centrality (truth selection applied)", HistType::kTH1F, {centAxis});

      histos.add("Correction/setSizes", "Sizes of sets", HistType::kTH1F, {{4, -0.5, 3.5}});
      auto hset = histos.get<TH1>(HIST("Correction/setSizes"));
      hset->GetXaxis()->SetBinLabel(1, "refClassIds");
      hset->GetXaxis()->SetBinLabel(2, "allowedMcIds");
      hset->GetXaxis()->SetBinLabel(3, "intersection");
      hset->GetXaxis()->SetBinLabel(4, "allowed-only");

      histos.add("Correction/hNEventsMCTruth", "hNEventsMCTruth", HistType::kTH1F, {AxisSpec{nSteps, 0.5, nSteps + 0.5, ""}});
      auto hstep = histos.get<TH1>(HIST("Correction/hNEventsMCTruth"));
      hstep->GetXaxis()->SetBinLabel(1, "All");
      hstep->GetXaxis()->SetBinLabel(2, "zvtx");
      hstep->GetXaxis()->SetBinLabel(3, "INEL>0");
      hstep->GetXaxis()->SetBinLabel(4, "Assoc with reco coll");

      histos.add("MCReco/hInvmass_Kstar_true", "MC-reco truth-tagged chK(892)", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});
      histos.add("MCReco/hInvmass_Kstar_bkg", "MC-reco residual background chK(892)", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});
    }

    ccdb->setURL(CCDBConfig.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in chK(892) Analysis Task";
    histos.print();
  }

  const int kCentFT0C = 1;
  const int kCentFT0M = 2;
  const float kInvalidCentrality = -999.f;

  template <typename CollisionType>
  float getCentrality(CollisionType const& collision)
  {
    if (cfgCentEst == kCentFT0C) {
      return collision.centFT0C();
    } else if (cfgCentEst == kCentFT0M) {
      return collision.centFT0M();
    } else {
      return kInvalidCentrality;
    }
  }

  // Track selection
  template <typename TrackType>
  bool trackCut(TrackType const& track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < TrackCuts.cfgMinPtcut)
      return false;
    if (std::abs(track.eta()) > TrackCuts.cfgMaxEtacut)
      return false;
    if (track.itsNCls() < TrackCuts.cfgITScluster)
      return false;
    if (track.tpcNClsFound() < TrackCuts.cfgTPCcluster)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < TrackCuts.cfgRatioTPCRowsOverFindableCls)
      return false;
    if (track.itsChi2NCl() >= TrackCuts.cfgITSChi2NCl)
      return false;
    if (track.tpcChi2NCl() >= TrackCuts.cfgTPCChi2NCl)
      return false;
    if (TrackCuts.cfgHasITS && !track.hasITS())
      return false;
    if (TrackCuts.cfgHasTPC && !track.hasTPC())
      return false;
    if (TrackCuts.cfgHasTOF && !track.hasTOF())
      return false;
    if (TrackCuts.cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (TrackCuts.cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    if (TrackCuts.cfgPVContributor && !track.isPVContributor())
      return false;
    if (TrackCuts.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (TrackCuts.cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (TrackCuts.cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (TrackCuts.cfgpTdepDCAxyCut) {
      // Tuned on the LHC22f anchored MC LHC23d1d on primary pions. 7 Sigmas of the resolution
      if (std::abs(track.dcaXY()) > (0.004 + (0.013 / track.pt())))
        return false;
    } else {
      if (std::abs(track.dcaXY()) > TrackCuts.cfgMaxbDCArToPVcut)
        return false;
    }
    if (std::abs(track.dcaZ()) > TrackCuts.cfgMaxbDCAzToPVcut)
      return false;
    return true;
  }

  // PID selection tools
  template <typename TrackType>
  bool selectionPIDPion(TrackType const& candidate)
  {
    if (std::abs(candidate.tpcNSigmaPi()) >= PIDCuts.cfgMaxTPCnSigmaPion)
      return false;
    if (PIDCuts.cfgTPConly)
      return true;
    if (candidate.pt() <= PIDCuts.cfgTOFMinPt)
      return true;

    if (candidate.hasTOF()) {
      const bool tofPIDPassed = std::abs(candidate.tofNSigmaPi()) < PIDCuts.cfgMaxTOFnSigmaPion;
      const bool combo = (PIDCuts.cfgNsigmaCutCombinedPion > 0) &&
                         (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() +
                            candidate.tofNSigmaPi() * candidate.tofNSigmaPi() <
                          PIDCuts.cfgNsigmaCutCombinedPion * PIDCuts.cfgNsigmaCutCombinedPion);
      return tofPIDPassed || combo;
    } else {
      return PIDCuts.cfgTOFVeto;
    }
  }

  template <typename CollisionType, typename K0sType>
  bool selectionK0s(CollisionType const& collision, K0sType const& candidate)
  {
    auto lDauDCA = std::fabs(candidate.dcaV0daughters());
    auto lDauPosDCAtoPV = std::fabs(candidate.dcapostopv());
    auto lDauNegDCAtoPV = std::fabs(candidate.dcanegtopv());
    auto lPt = candidate.pt();
    auto lRapidity = candidate.yK0Short();
    auto lRadius = candidate.v0radius();
    auto lDCAtoPV = std::fabs(candidate.dcav0topv());
    auto lCPA = candidate.v0cosPA();
    auto lPropTauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassK0Short;
    auto lMk0s = candidate.mK0Short();
    auto lMLambda = candidate.mLambda();
    auto lMALambda = candidate.mAntiLambda();

    auto checkCommonCuts = [&]() {
      if (std::fabs(lDauDCA) > SecondaryCuts.cfgSecondaryDauDCAMax)
        return false;
      if (std::fabs(lDauPosDCAtoPV) < SecondaryCuts.cfgSecondaryDauPosDCAtoPVMin)
        return false;
      if (std::fabs(lDauNegDCAtoPV) < SecondaryCuts.cfgSecondaryDauNegDCAtoPVMin)
        return false;
      if (lPt < SecondaryCuts.cfgSecondaryPtMin)
        return false;
      if (std::fabs(lRapidity) > SecondaryCuts.cfgSecondaryRapidityMax)
        return false;
      if (lRadius < SecondaryCuts.cfgSecondaryRadiusMin || lRadius > SecondaryCuts.cfgSecondaryRadiusMax)
        return false;
      if (std::fabs(lDCAtoPV) > SecondaryCuts.cfgSecondaryDCAtoPVMax)
        return false;
      if (lCPA < SecondaryCuts.cfgSecondaryCosPAMin)
        return false;
      if (lPropTauK0s > SecondaryCuts.cfgSecondaryProperLifetimeMax)
        return false;
      if (candidate.qtarm() < SecondaryCuts.cfgSecondaryparamArmenterosCut * std::fabs(candidate.alpha()))
        return false;
      if (std::fabs(lMk0s - MassK0Short) > SecondaryCuts.cfgSecondaryMassWindow)
        return false;
      if (SecondaryCuts.cfgSecondaryCrossMassHypothesisCut &&
          ((std::fabs(lMLambda - MassLambda0) < SecondaryCuts.cfgSecondaryCrossMassCutWindow) || (std::fabs(lMALambda - MassLambda0Bar) < SecondaryCuts.cfgSecondaryCrossMassCutWindow)))
        return false;
      return true;
    };

    if (SecondaryCuts.cfgReturnFlag) { // For cut study
      bool returnFlag = true;
      histos.fill(HIST("QA/K0sCutCheck"), 0);
      if (lDauDCA > SecondaryCuts.cfgSecondaryDauDCAMax) {
        histos.fill(HIST("QA/K0sCutCheck"), 1);
        returnFlag = false;
      }
      if (lDauPosDCAtoPV < SecondaryCuts.cfgSecondaryDauPosDCAtoPVMin) {
        histos.fill(HIST("QA/K0sCutCheck"), 2);
        returnFlag = false;
      }
      if (lDauNegDCAtoPV < SecondaryCuts.cfgSecondaryDauNegDCAtoPVMin) {
        histos.fill(HIST("QA/K0sCutCheck"), 3);
        returnFlag = false;
      }
      if (lPt < SecondaryCuts.cfgSecondaryPtMin) {
        histos.fill(HIST("QA/K0sCutCheck"), 4);
        returnFlag = false;
      }
      if (std::fabs(lRapidity) > SecondaryCuts.cfgSecondaryRapidityMax) {
        histos.fill(HIST("QA/K0sCutCheck"), 5);
        returnFlag = false;
      }
      if (lRadius < SecondaryCuts.cfgSecondaryRadiusMin || lRadius > SecondaryCuts.cfgSecondaryRadiusMax) {
        histos.fill(HIST("QA/K0sCutCheck"), 6);
        returnFlag = false;
      }
      if (lDCAtoPV > SecondaryCuts.cfgSecondaryDCAtoPVMax) {
        histos.fill(HIST("QA/K0sCutCheck"), 7);
        returnFlag = false;
      }
      if (lCPA < SecondaryCuts.cfgSecondaryCosPAMin) {
        histos.fill(HIST("QA/K0sCutCheck"), 8);
        returnFlag = false;
      }
      if (lPropTauK0s > SecondaryCuts.cfgSecondaryProperLifetimeMax) {
        histos.fill(HIST("QA/K0sCutCheck"), 9);
        returnFlag = false;
      }
      if (candidate.qtarm() < SecondaryCuts.cfgSecondaryparamArmenterosCut * std::fabs(candidate.alpha())) {
        histos.fill(HIST("QA/K0sCutCheck"), 10);
        returnFlag = false;
      }
      if (std::fabs(lMk0s - MassK0Short) > SecondaryCuts.cfgSecondaryMassWindow) {
        histos.fill(HIST("QA/K0sCutCheck"), 11);
        returnFlag = false;
      }
      if (SecondaryCuts.cfgSecondaryCrossMassHypothesisCut &&
          ((std::fabs(lMLambda - MassLambda0) < SecondaryCuts.cfgSecondaryCrossMassCutWindow) || (std::fabs(lMALambda - MassLambda0Bar) < SecondaryCuts.cfgSecondaryCrossMassCutWindow))) {
        histos.fill(HIST("QA/K0sCutCheck"), 12);
        returnFlag = false;
      }
      return returnFlag;
    } else { // normal usage
      if (SecondaryCuts.cfgSecondaryRequire) {
        return checkCommonCuts();
      } else {
        return std::fabs(lMk0s - MassK0Short) <= SecondaryCuts.cfgSecondaryMassWindow; // always apply mass window cut
      }
    }
  } // selectionK0s

  template <typename TrackTemplate, typename V0Template>
  bool isTrueKstar(const TrackTemplate& bTrack, const V0Template& K0scand)
  {
    if (std::abs(bTrack.PDGCode()) != kPiPlus) // Are you pion?
      return false;
    if (std::abs(K0scand.PDGCode()) != kPDGK0s) // Are you K0s?
      return false;

    auto motherbTrack = bTrack.template mothers_as<aod::McParticles>();
    auto motherkV0 = K0scand.template mothers_as<aod::McParticles>();

    // Check bTrack first
    if (std::abs(motherbTrack.pdgCode()) != kKstarPlus) // Are you charged Kstar's daughter?
      return false;                                     // Apply first since it's more restrictive

    if (std::abs(motherkV0.pdgCode()) != kPDGK0s) // Is it K0s?
      return false;
    // Check if K0s's mother is K0 (311)
    auto motherK0 = motherkV0.template mothers_as<aod::McParticles>();
    if (std::abs(motherK0.pdgCode()) != kPDGK0)
      return false;

    // Check if K0's mother is Kstar (323)
    auto motherKstar = motherK0.template mothers_as<aod::McParticles>();
    if (std::abs(motherKstar.pdgCode()) != kKstarPlus)
      return false;

    // Check if bTrack and K0 have the same mother (global index)
    if (motherbTrack.globalIndex() != motherK0.globalIndex())
      return false;

    return true;
  }

  std::unordered_set<int64_t> allowedMcIds;
  std::unordered_map<int64_t, float> centTruthByAllowed;
  std::unordered_set<int64_t> refClassIds;
  std::unordered_map<int64_t, float> refCentByMcId;

  template <typename RecoEventsT>
  void buildAllowedMcIds(RecoEventsT const& events)
  {
    allowedMcIds.clear();
    centTruthByAllowed.clear();

    for (const auto& coll : events) {
      // lCentrality = getCentrality(coll);

      if (!coll.has_mcCollision())
        continue;

      const auto mcid = coll.mcCollisionId();
      const auto mccoll = coll.template mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();

      const float lCentrality = mccoll.centFT0M();

      if (doprocessMC) {
        histos.fill(HIST("QACent_woCut"), lCentrality);
        histos.fill(HIST("QAvtxz_woCut"), coll.posZ());
      }

      if (!colCuts.isSelected(coll))
        continue;
      if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(coll))
        continue;
      if (!coll.isInelGt0())
        continue;

      if (doprocessMC) {
        histos.fill(HIST("QACent_woCentCut"), lCentrality);
        histos.fill(HIST("QAvtxz_wVtxzCut"), coll.posZ());
      }

      if (lCentrality < EventCuts.cfgEventCentralityMin || lCentrality > EventCuts.cfgEventCentralityMax)
        continue;

      if (doprocessMC) {
        histos.fill(HIST("QACent_wCentCut"), lCentrality);
      }
      allowedMcIds.insert(mcid);
      centTruthByAllowed.emplace(mcid, lCentrality);
    }
  }

  template <typename McCollsT, typename McPartsT>
  void buildReferenceMcIds(McCollsT const& mccolls, McPartsT const& mcparts)
  {
    refClassIds.clear();
    refCentByMcId.clear();

    for (const auto& coll : mccolls) {
      bool pass = true;

      if (cfgTruthIncludeZvtx && std::abs(coll.posZ()) >= EventCuts.cfgEvtZvtx)
        pass = false;

      if (pass && cfgTruthUseInelGt0) {
        auto partsThisMc = mcparts.sliceBy(perMCCollision, coll.globalIndex());
        if (!pwglf::isINELgtNmc(partsThisMc, 0, pdg))
          pass = false;
      }

      if (!pass)
        continue;

      const auto mcid = coll.globalIndex();
      refClassIds.insert(mcid);
      const float lCentrality = coll.centFT0M();
      refCentByMcId.emplace(mcid, lCentrality);
    }
  }

  void effK0sProcessGen(MCTrueTrackCandidates const& mcparts)
  {
    for (const auto& part : mcparts) {
      if (!part.has_mcCollision())
        continue;
      if (part.pdgCode() != kPDGK0s)
        continue;
      if (!part.isPhysicalPrimary())
        continue;
      if (std::abs(part.y()) > SecondaryCuts.cfgSecondaryRapidityMax)
        continue;

      const auto mcid = part.mcCollisionId();
      if (allowedMcIds.count(mcid) == 0)
        continue;

      auto iter = centTruthByAllowed.find(mcid);
      if (iter == centTruthByAllowed.end())
        continue;

      const float lCentrality = iter->second;

      histos.fill(HIST("EffK0s/genK0s"), part.pt(), lCentrality);
    }
  }

  void effK0sProcessReco(MCV0Candidates const& v0s)
  {
    for (const auto& v0 : v0s) {
      auto coll = v0.template collision_as<MCEventCandidates>();

      if (!coll.has_mcCollision())
        continue;

      const auto mcid = coll.mcCollisionId();

      if (allowedMcIds.count(mcid) == 0)
        continue;

      const auto mccoll = coll.template mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      const float lCentrality = mccoll.centFT0M();

      const double ptreco = v0.pt();
      const double yreco = v0.yK0Short();

      double ptgen = -1, ygen = 0;
      if (!getTruthK0sAndGenKinematics(v0, ptgen, ygen))
        continue;
      if (std::abs(yreco) > SecondaryCuts.cfgSecondaryRapidityMax)
        continue;

      if (!SecondaryCuts.cfgByPassDauPIDSelection) {
        auto posDauTrack = v0.template posTrack_as<MCTrackCandidates>();
        auto negDauTrack = v0.template negTrack_as<MCTrackCandidates>();
        if (!selectionPIDPion(posDauTrack))
          continue;
        if (!selectionPIDPion(negDauTrack))
          continue;
      }
      if (!selectionK0s(coll, v0))
        continue;

      histos.fill(HIST("EffK0s/recoK0s"), ptreco, lCentrality);
    }
  } // effK0sProcessReco

  template <typename V0T, typename TrkT>
  bool matchRecoToTruthKstar(V0T const& v0, TrkT const& trk, double& ptgen, double& ygen)
  {
    if (!v0.has_mcParticle() || !trk.has_mcParticle())
      return false;

    auto mcK0s = v0.template mcParticle_as<MCTrueTrackCandidates>();
    auto mcPi = trk.template mcParticle_as<MCTrueTrackCandidates>();

    if (std::abs(mcK0s.pdgCode()) != kPDGK0s)
      return false;
    if (std::abs(mcPi.pdgCode()) != kPiPlus)
      return false;

    MCTrueTrackCandidates::iterator kstarFromPi;
    bool havePiKstar = false;
    for (const auto& m1 : mcPi.template mothers_as<MCTrueTrackCandidates>()) {
      if (std::abs(m1.pdgCode()) == kKstarPlus) {
        kstarFromPi = m1;
        havePiKstar = true;
        break;
      }
    }
    if (!havePiKstar) {
      return false;
    }

    bool shareSameKstar = false;
    for (const auto& m1 : mcK0s.template mothers_as<MCTrueTrackCandidates>()) {
      if (std::abs(m1.pdgCode()) == kPDGK0) {
        for (const auto& m2 : m1.template mothers_as<MCTrueTrackCandidates>()) {
          if (m2.globalIndex() == kstarFromPi.globalIndex()) {
            shareSameKstar = true;
            break;
          }
        }
        if (shareSameKstar)
          break;
      }
    }
    if (!shareSameKstar) {
      return false;
    }

    ptgen = kstarFromPi.pt();
    ygen = kstarFromPi.y();

    return true;
  } // matchRecoToTruthKstar

  void effKstarProcessGen(MCTrueTrackCandidates const& mcparts)
  {
    for (const auto& part : mcparts) {
      if (!part.has_mcCollision())
        continue;
      if (std::abs(part.pdgCode()) != kKstarPlus)
        continue;
      if (std::abs(part.y()) > KstarCuts.cfgKstarMaxRap)
        continue;

      const auto mcid = part.mcCollisionId();
      if (allowedMcIds.count(mcid) == 0)
        continue;

      auto iter = centTruthByAllowed.find(mcid);
      if (iter == centTruthByAllowed.end())
        continue;

      const float lCentrality = iter->second;

      histos.fill(HIST("EffKstar/genKstar"), part.pt(), lCentrality);
    }
  } // effKstarProcessGen

  template <typename V0RangeT, typename TrkRangeT>
  void effKstarProcessReco(V0RangeT const& v0s, TrkRangeT const& tracks)
  {
    for (const auto& v0 : v0s) {
      auto coll = v0.template collision_as<MCEventCandidates>();

      if (!coll.has_mcCollision())
        continue;

      const auto mcid = coll.mcCollisionId();

      if (allowedMcIds.count(mcid) == 0)
        continue;

      const auto mccoll = coll.template mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      const float lCentrality = mccoll.centFT0M();

      if (!SecondaryCuts.cfgByPassDauPIDSelection) {
        auto posDauTrack = v0.template posTrack_as<MCTrackCandidates>();
        auto negDauTrack = v0.template negTrack_as<MCTrackCandidates>();
        if (!selectionPIDPion(posDauTrack))
          continue;
        if (!selectionPIDPion(negDauTrack))
          continue;
      }
      if (!selectionK0s(coll, v0))
        continue;

      auto trks = tracks.sliceBy(perCollision, v0.collisionId());
      for (const auto& bTrack : trks) {
        if (bTrack.collisionId() != v0.collisionId())
          continue;
        if (!trackCut(bTrack))
          continue;
        if (!selectionPIDPion(bTrack))
          continue;

        LorentzVectorSetXYZM lResoSecondary, lDecayDaughter_bach, lResoKstar, lDaughterRot, lResonanceRot;

        lResoSecondary = LorentzVectorSetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
        lDecayDaughter_bach = LorentzVectorSetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), MassPionCharged);
        lResoKstar = lResoSecondary + lDecayDaughter_bach;

        const double ptreco = lResoKstar.Pt();
        const double yreco = lResoKstar.Rapidity();

        if (std::abs(yreco) > KstarCuts.cfgKstarMaxRap)
          continue;

        double ptgen = 0, ygen = 0;
        const bool isTrue = matchRecoToTruthKstar(v0, bTrack, ptgen, ygen);

        if (isTrue) {

          histos.fill(HIST("EffKstar/recoKstar"), ptreco, lCentrality);
          histos.fill(HIST("MCReco/hInvmass_Kstar_true"), lCentrality, ptreco, lResoKstar.M());

        } else {
          histos.fill(HIST("MCReco/hInvmass_Kstar_bkg"), lCentrality, ptreco, lResoKstar.M());
        }
      }
    }
  } // effKstarProcessReco

  void fillSigLossNum(MCTrueTrackCandidates const& mcparts)
  {
    for (auto const& part : mcparts) {
      if (!part.has_mcCollision())
        continue;
      if (std::abs(part.pdgCode()) != kKstarPlus)
        continue;
      if (std::abs(part.y()) > KstarCuts.cfgKstarMaxRap)
        continue;

      const auto mcid = part.mcCollisionId();
      if (allowedMcIds.count(mcid) == 0)
        continue;

      auto iter = centTruthByAllowed.find(mcid);
      if (iter == centTruthByAllowed.end())
        continue;

      const float lCentrality = iter->second;

      histos.fill(HIST("Correction/sigLoss_num"), part.pt(), lCentrality);
    }
  } // fillSigLossNum

  void fillSigLossDen(MCTrueTrackCandidates const& mcparts)
  {
    for (auto const& part : mcparts) {
      if (!part.has_mcCollision())
        continue;
      if (std::abs(part.pdgCode()) != kKstarPlus)
        continue;
      if (std::abs(part.y()) > KstarCuts.cfgKstarMaxRap)
        continue;

      const auto mcid = part.mcCollisionId();
      if (refClassIds.count(mcid) == 0)
        continue;

      auto iter = refCentByMcId.find(mcid);
      if (iter == refCentByMcId.end())
        continue;

      const float lCentrality = iter->second;

      histos.fill(HIST("Correction/sigLoss_den"), part.pt(), lCentrality);
    }
  } // fillSigLossDen

  int count = 0;

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType, typename TracksTypeK0s>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksTypeK0s& dTracks2)
  {
    using TrackTarget = std::decay_t<TracksType>;

    histos.fill(HIST("QA/before/CentDist"), lCentrality);

    LorentzVectorSetXYZM lDecayDaughter1, lDecayDaughter2, lResoSecondary, lDecayDaughter_bach, lResoKstar, lDaughterRot, lResonanceRot;
    std::vector<int> trackIndicies = {};
    std::vector<int> k0sIndicies = {};

    for (const auto& bTrack : dTracks1) {
      auto trkbpt = bTrack.pt();
      auto istrkbhasTOF = bTrack.hasTOF();
      auto trkbNSigmaPiTPC = bTrack.tpcNSigmaPi();
      auto trkbNSigmaPiTOF = (istrkbhasTOF) ? bTrack.tofNSigmaPi() : -999.;

      if constexpr (!IsMix) {
        if (cfgFillQAPlots) {
          // Bachelor pion QA plots
          histos.fill(HIST("QA/before/trkbpionTPCPID"), trkbpt, trkbNSigmaPiTPC);
          if (istrkbhasTOF) {
            histos.fill(HIST("QA/before/trkbpionTOFPID"), trkbpt, trkbNSigmaPiTOF);
            histos.fill(HIST("QA/before/trkbpionTPCTOFPID"), trkbNSigmaPiTPC, trkbNSigmaPiTOF);
          }
          histos.fill(HIST("QA/before/trkbpionpT"), trkbpt);
          histos.fill(HIST("QA/before/trkbpionDCAxy"), bTrack.dcaXY());
          histos.fill(HIST("QA/before/trkbpionDCAz"), bTrack.dcaZ());
        }
      }

      if (!trackCut(bTrack))
        continue;
      if (!selectionPIDPion(bTrack))
        continue;

      if constexpr (!IsMix) {
        if (cfgFillQAPlots) {
          // Bachelor pion QA plots after applying cuts
          histos.fill(HIST("QA/after/trkbpionTPCPID"), trkbpt, trkbNSigmaPiTPC);
          if (istrkbhasTOF) {
            histos.fill(HIST("QA/after/trkbpionTOFPID"), trkbpt, trkbNSigmaPiTOF);
            histos.fill(HIST("QA/after/trkbpionTPCTOFPID"), trkbNSigmaPiTPC, trkbNSigmaPiTOF);
          }
          histos.fill(HIST("QA/after/trkbpionpT"), trkbpt);
          histos.fill(HIST("QA/after/trkbpionDCAxy"), bTrack.dcaXY());
          histos.fill(HIST("QA/after/trkbpionDCAz"), bTrack.dcaZ());
        }
      }
      trackIndicies.push_back(bTrack.index());
    }

    for (const auto& k0sCand : dTracks2) {

      auto posDauTrack = k0sCand.template posTrack_as<TrackTarget>();
      auto negDauTrack = k0sCand.template negTrack_as<TrackTarget>();

      /// Daughters
      // Positve pion
      auto trkppt = posDauTrack.pt();
      auto istrkphasTOF = posDauTrack.hasTOF();
      auto trkpNSigmaPiTPC = posDauTrack.tpcNSigmaPi();
      auto trkpNSigmaPiTOF = (istrkphasTOF) ? posDauTrack.tofNSigmaPi() : -999.;
      // Negative pion
      auto trknpt = negDauTrack.pt();
      auto istrknhasTOF = negDauTrack.hasTOF();
      auto trknNSigmaPiTPC = negDauTrack.tpcNSigmaPi();
      auto trknNSigmaPiTOF = (istrknhasTOF) ? negDauTrack.tofNSigmaPi() : -999.;

      /// K0s
      auto trkkDauDCA = k0sCand.dcaV0daughters();
      auto trkky = k0sCand.yK0Short();
      auto trkkDCAtoPV = k0sCand.dcav0topv();
      auto trkkCPA = k0sCand.v0cosPA();
      auto trkkPropTau = k0sCand.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassK0Short;
      auto trkkMass = k0sCand.mK0Short();
      auto trkkDauDCAPostoPV = k0sCand.dcapostopv();
      auto trkkDauDCANegtoPV = k0sCand.dcanegtopv();
      auto trkkpt = k0sCand.pt();
      auto trkkRadius = k0sCand.v0radius();

      lResoSecondary = LorentzVectorSetXYZM(k0sCand.px(), k0sCand.py(), k0sCand.pz(), trkkMass);

      if constexpr (!IsMix) {
        if (cfgFillQAPlots) {
          // Seconddary QA plots
          histos.fill(HIST("QA/before/trkppionTPCPID"), trkppt, trkpNSigmaPiTPC);
          if (istrkphasTOF) {
            histos.fill(HIST("QA/before/trkppionTOFPID"), trkppt, trkpNSigmaPiTOF);
            histos.fill(HIST("QA/before/trkppionTPCTOFPID"), trkpNSigmaPiTPC, trkpNSigmaPiTOF);
          }
          histos.fill(HIST("QA/before/trkppionpT"), trkppt);
          histos.fill(HIST("QA/before/trkppionDCAxy"), posDauTrack.dcaXY());
          histos.fill(HIST("QA/before/trkppionDCAz"), posDauTrack.dcaZ());

          histos.fill(HIST("QA/before/trknpionTPCPID"), trknpt, trknNSigmaPiTPC);
          if (istrknhasTOF) {
            histos.fill(HIST("QA/before/trknpionTOFPID"), trknpt, trknNSigmaPiTOF);
            histos.fill(HIST("QA/before/trknpionTPCTOFPID"), trknNSigmaPiTPC, trknNSigmaPiTOF);
          }
          histos.fill(HIST("QA/before/trknpionpT"), trknpt);
          histos.fill(HIST("QA/before/trknpionDCAxy"), negDauTrack.dcaXY());
          histos.fill(HIST("QA/before/trknpionDCAz"), negDauTrack.dcaZ());

          histos.fill(HIST("QA/before/hDauDCASecondary"), trkkDauDCA);
          histos.fill(HIST("QA/before/hDauPosDCAtoPVSecondary"), trkkDauDCAPostoPV);
          histos.fill(HIST("QA/before/hDauNegDCAtoPVSecondary"), trkkDauDCANegtoPV);

          histos.fill(HIST("QA/before/hpT_Secondary"), trkkpt);
          histos.fill(HIST("QA/before/hy_Secondary"), trkky);
          histos.fill(HIST("QA/before/hRadiusSecondary"), trkkRadius);
          histos.fill(HIST("QA/before/hDCAtoPVSecondary"), trkkDCAtoPV);
          histos.fill(HIST("QA/before/hCPASecondary"), trkkCPA);
          histos.fill(HIST("QA/before/hPropTauSecondary"), trkkPropTau);
          histos.fill(HIST("QA/before/hInvmassSecondary"), trkkMass);
        }
      }

      if (!SecondaryCuts.cfgByPassDauPIDSelection && !selectionPIDPion(posDauTrack))
        continue;
      if (!SecondaryCuts.cfgByPassDauPIDSelection && !selectionPIDPion(negDauTrack))
        continue;
      if (!selectionK0s(collision, k0sCand))
        continue;

      if constexpr (!IsMix) {
        if (cfgFillQAPlots) {
          // Seconddary QA plots after applying cuts

          histos.fill(HIST("QA/after/trkppionTPCPID"), trkppt, trkpNSigmaPiTPC);
          if (istrkphasTOF) {
            histos.fill(HIST("QA/after/trkppionTOFPID"), trkppt, trkpNSigmaPiTOF);
            histos.fill(HIST("QA/after/trkppionTPCTOFPID"), trkpNSigmaPiTPC, trkpNSigmaPiTOF);
          }
          histos.fill(HIST("QA/after/trkppionpT"), trkppt);
          histos.fill(HIST("QA/after/trkppionDCAxy"), posDauTrack.dcaXY());
          histos.fill(HIST("QA/after/trkppionDCAz"), posDauTrack.dcaZ());

          histos.fill(HIST("QA/after/trknpionTPCPID"), trknpt, trknNSigmaPiTPC);
          if (istrknhasTOF) {
            histos.fill(HIST("QA/after/trknpionTOFPID"), trknpt, trknNSigmaPiTOF);
            histos.fill(HIST("QA/after/trknpionTPCTOFPID"), trknNSigmaPiTPC, trknNSigmaPiTOF);
          }
          histos.fill(HIST("QA/after/trknpionpT"), trknpt);
          histos.fill(HIST("QA/after/trknpionDCAxy"), negDauTrack.dcaXY());
          histos.fill(HIST("QA/after/trknpionDCAz"), negDauTrack.dcaZ());

          histos.fill(HIST("QA/after/hDauDCASecondary"), trkkDauDCA);
          histos.fill(HIST("QA/after/hDauPosDCAtoPVSecondary"), trkkDauDCAPostoPV);
          histos.fill(HIST("QA/after/hDauNegDCAtoPVSecondary"), trkkDauDCANegtoPV);

          histos.fill(HIST("QA/after/hpT_Secondary"), trkkpt);
          histos.fill(HIST("QA/after/hy_Secondary"), trkky);
          histos.fill(HIST("QA/after/hRadiusSecondary"), trkkRadius);
          histos.fill(HIST("QA/after/hDCAtoPVSecondary"), trkkDCAtoPV);
          histos.fill(HIST("QA/after/hCPASecondary"), trkkCPA);
          histos.fill(HIST("QA/after/hPropTauSecondary"), trkkPropTau);
          histos.fill(HIST("QA/after/hInvmassSecondary"), trkkMass);
        }
        histos.fill(HIST("hInvmass_K0s"), lCentrality, lResoSecondary.Pt(), lResoSecondary.M());
      }
      k0sIndicies.push_back(k0sCand.index());
    }

    for (const auto& trackIndex : trackIndicies) {
      auto bTrack = dTracks1.rawIteratorAt(trackIndex);
      for (const auto& k0sIndex : k0sIndicies) {
        auto k0sCand = dTracks2.rawIteratorAt(k0sIndex);

        auto trkkMass = k0sCand.mK0Short();

        lDecayDaughter_bach = LorentzVectorSetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), MassPionCharged);
        lResoSecondary = LorentzVectorSetXYZM(k0sCand.px(), k0sCand.py(), k0sCand.pz(), trkkMass);
        lResoKstar = lResoSecondary + lDecayDaughter_bach;

        // QA plots
        if constexpr (!IsMix) {
          if (cfgFillQAPlots) {
            histos.fill(HIST("QA/before/KstarRapidity"), lResoKstar.Rapidity());
            histos.fill(HIST("QA/before/kstarinvmass"), lResoKstar.M());
          }
        }

        if (lResoKstar.Rapidity() > KstarCuts.cfgKstarMaxRap || lResoKstar.Rapidity() < KstarCuts.cfgKstarMinRap)
          continue;

        if constexpr (!IsMix) {
          unsigned int typeKstar = bTrack.sign() > 0 ? BinType::kKstarP : BinType::kKstarN;
          if (cfgFillQAPlots) {

            histos.fill(HIST("QA/after/KstarRapidity"), lResoKstar.Rapidity());
            histos.fill(HIST("QA/after/kstarinvmass"), lResoKstar.M());
          }
          histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResoKstar.Pt(), lResoKstar.M());

          if (BkgEstimationConfig.cfgFillRotBkg) {
            for (int i = 0; i < BkgEstimationConfig.cfgNrotBkg; i++) {
              auto lRotAngle = BkgEstimationConfig.cfgMinRot + i * ((BkgEstimationConfig.cfgMaxRot - BkgEstimationConfig.cfgMinRot) / (BkgEstimationConfig.cfgNrotBkg - 1));
              if (cfgFillQAPlots) {
                histos.fill(HIST("QA/RotBkg/hRotBkg"), lRotAngle);
              }
              if (BkgEstimationConfig.cfgRotPion) {
                lDaughterRot = lDecayDaughter_bach;
                // lDaughterRot.RotateZ(lRotAngle);
                ROOT::Math::RotationZ rot(lRotAngle);
                auto p3 = rot * lDaughterRot.Vect();
                lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                lResonanceRot = lDaughterRot + lResoSecondary;
              } else {
                lDaughterRot = lResoSecondary;
                // lDaughterRot.RotateZ(lRotAngle);
                ROOT::Math::RotationZ rot(lRotAngle);
                auto p3 = rot * lDaughterRot.Vect();
                lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                lResonanceRot = lDecayDaughter_bach + lDaughterRot;
              }
              typeKstar = bTrack.sign() > 0 ? BinType::kKstarP_Rot : BinType::kKstarN_Rot;
              histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResonanceRot.Pt(), lResonanceRot.M());
            }
          }
        } // IsMix
      } // K0scand
    } // bTrack

    count++;

  } // fillHistograms

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks,
                   V0Candidates const& v0s,
                   aod::BCsWithTimestamps const&)
  {
    if (!colCuts.isSelected(collision)) // Default event selection
      return;
    if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(collision)) {
      return;
    }
    lCentrality = getCentrality(collision);
    if (lCentrality < EventCuts.cfgEventCentralityMin || lCentrality > EventCuts.cfgEventCentralityMax)
      return;
    if (!collision.isInelGt0())
      return;
    colCuts.fillQA(collision);

    fillHistograms<false, false>(collision, tracks, v0s);
  }
  PROCESS_SWITCH(Chk892pp, processData, "Process Event for data without Partitioning", false);

  void processMC(MCTrueTrackCandidates const& mcpart,
                 MCTrackCandidates const& tracks,
                 MCV0Candidates const& v0s,
                 MCEventCandidates const& events,
                 soa::Join<MCTrueEventCandidates, aod::McCentFT0Ms> const& mccolls)
  {
    buildAllowedMcIds(events);
    buildReferenceMcIds(mccolls, mcpart);
    effK0sProcessGen(mcpart);
    effK0sProcessReco(v0s);
    effKstarProcessGen(mcpart);
    effKstarProcessReco(v0s, tracks);
    fillSigLossNum(mcpart);
    fillSigLossDen(mcpart);

    for (const auto& mcid : refClassIds) {
      histos.fill(HIST("Correction/EF_den"), refCentByMcId[mcid]);
    }
    for (const auto& mcid : allowedMcIds) {
      auto iter = centTruthByAllowed.find(mcid);
      if (iter == centTruthByAllowed.end())
        continue;

      const float lCentrality = iter->second;
      histos.fill(HIST("Correction/EF_num"), lCentrality);
    }

    size_t nIntersect = 0;
    for (const auto& mcid : allowedMcIds)
      if (refClassIds.count(mcid))
        nIntersect++;
    histos.fill(HIST("Correction/setSizes"), 0.0, refClassIds.size());
    histos.fill(HIST("Correction/setSizes"), 1.0, allowedMcIds.size());
    histos.fill(HIST("Correction/setSizes"), 2.0, nIntersect);
    histos.fill(HIST("Correction/setSizes"), 3.0, allowedMcIds.size() - nIntersect);

    for (const auto& mcc : mccolls) {
      histos.fill(HIST("Correction/MCTruthCent_all"), mcc.centFT0M());
    }

    for (const auto& mcid : refClassIds) {
      auto iter = refCentByMcId.find(mcid);
      if (iter == refCentByMcId.end())
        continue;
      lCentrality = iter->second;
      histos.fill(HIST("Correction/MCTruthCent_cut"), lCentrality);
    }

    for (auto const& mcc : mccolls) {
      const auto mcid = mcc.globalIndex();

      histos.fill(HIST("Correction/hNEventsMCTruth"), 1.0);

      bool passZvtx = true;
      if (cfgTruthIncludeZvtx && std::abs(mcc.posZ()) > EventCuts.cfgEvtZvtx) {
        passZvtx = false;
      }
      if (passZvtx) {
        histos.fill(HIST("Correction/hNEventsMCTruth"), 2.0);

        auto partsThisMc = mcpart.sliceBy(perMCCollision, mcid);
        if (pwglf::isINELgtNmc(partsThisMc, 0, pdg)) {
          histos.fill(HIST("Correction/hNEventsMCTruth"), 3.0);
        }
      }
      if (allowedMcIds.count(mcid)) {
        histos.fill(HIST("Correction/hNEventsMCTruth"), 4.0);
      }
    }
  }
  PROCESS_SWITCH(Chk892pp, processMC, "Process Event for MC", true);

  void processMCQA(MCEventCandidates::iterator const& collision,
                   MCTrackCandidates const& tracks,
                   MCV0Candidates const& v0s,
                   soa::Join<MCTrueEventCandidates, aod::McCentFT0Ms> const& mccolls,
                   aod::BCsWithTimestamps const&)
  {
    if (!colCuts.isSelected(collision))
      return;
    if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(collision))
      return;
    if (!collision.isInelGt0())
      return;

    if (!collision.has_mcCollision())
      return;

    auto id = collision.mcCollisionId();

    auto mccoll = mccolls.iteratorAt(id);
    const float lCentrality = mccoll.centFT0M();

    if (lCentrality < EventCuts.cfgEventCentralityMin || lCentrality > EventCuts.cfgEventCentralityMax)
      return;
    colCuts.fillQA(collision);

    fillHistograms<true, false>(collision, tracks, v0s);
  }
  PROCESS_SWITCH(Chk892pp, processMCQA, "Process Event for MC and fill QA plots", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Chk892pp>(cfgc)};
}
