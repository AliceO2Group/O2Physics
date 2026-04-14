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

/// \file k1analysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>

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
using namespace o2::constants::math;
using namespace o2::aod::rctsel;

struct K1analysis {
  enum BinType : unsigned int {
    kK1P = 0,
    kK1A,
    kK1P_Like,
    kK1A_Like,
    kK1P_Rot,
    kK1A_Rot,
    kTYEnd
  };

  enum EvtStep {
    kAll = 0,
    kZvtx,
    kINELgt0,
    kAssocReco,
    kNSteps
  };

  enum class K1MassRegion : uint8_t {
    Outside = 0,
    Signal,
    SBLeft,
    SBRight
  };

  const int nSteps = static_cast<int>(EvtStep::kNSteps);

  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  // Preslice<aod::V0s> perCollisionV0 = aod::v0data::collisionId;
  Preslice<aod::V0Datas> perCollisionV0 = aod::v0data::collisionId;
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
    ConfigurableAxis cfgBinsPt{"cfgBinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsPtQA{"cfgBinsPtQA", {VARIABLE_WIDTH, 0.0, 0.3, 0.6, 1.2, 1.8, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0, 10.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsCent{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 110.0}, "Binning of the centrality axis"};
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

  Configurable<bool> cfgFillQAPlots{"cfgFillQAPlots", true, "Fill QA plots"};
  Configurable<int> cfgCentEst{"cfgCentEst", 2, "Centrality estimator, 1: FT0C, 2: FT0M"};

  /// PID Selections, pion
  struct : ConfigurableGroup {
    Configurable<bool> cfgTPConly{"cfgTPConly", true, "Use only TPC for PID"};                                      // bool
    Configurable<float> cfgMaxTPCnSigmaPion{"cfgMaxTPCnSigmaPion", 5.0, "TPC nSigma cut for Pion"};                 // TPC
    Configurable<float> cfgMaxTOFnSigmaPion{"cfgMaxTOFnSigmaPion", 5.0, "TOF nSigma cut for Pion"};                 // TOF
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
    Configurable<bool> cfgpTdepDCAzCut{"cfgpTdepDCAzCut", true, "pT-dependent DCAz cut"};
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
    Configurable<float> cfgKstarMinPt{"cfgKstarMinPt", 0.0, "Kstar minimum pT"};
    Configurable<float> cfgKstarMaxRap{"cfgKstarMaxRap", 0.5, "Kstar maximum rapidity"};
    Configurable<float> cfgKstarMinRap{"cfgKstarMinRap", -0.5, "Kstar minimum rapidity"};
    Configurable<float> cfgKstarMassWindow{"cfgKstarMassWindow", 0.05, "Kstar inv mass selection window"};
  } KstarCuts;

  // K1 selection
  struct : ConfigurableGroup {
    Configurable<double> cfgK1MinPt{"cfgK1MinPt", 0.0, "K1 minimum pT"};
    Configurable<double> cfgK1MaxRap{"cfgK1MaxRap", 0.5, "K1 maximum rapidity"};
    Configurable<double> cfgK1MinRap{"cfgK1MinRap", -0.5, "K1 minimum rapidity"};
  } K1Cuts;

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
  Configurable<bool> cfgHasPair{"cfgHasPair", true, "Check the existence of pairs"};
  Configurable<float> cfgPiPiMinPt{"cfgPiPiMinPt", 0.5, "bachelor pion and secondary pion minimum pT cut"};

  float lCentrality;

  // PDG code
  int kPDGK0s = kK0Short;
  int kPDGK0 = kK0;
  int kPDGKstarPlus = o2::constants::physics::Pdg::kKPlusStar892;
  int kPDGPiPlus = kPiPlus;
  int kPDGK10 = 10313;
  double fMaxPosPV = 1e-2;

  void init(o2::framework::InitContext&)
  {
    lCentrality = -999;

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
    AxisSpec ptAxis = {AxisConfig.cfgBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {AxisConfig.cfgBinsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec radiusAxis = {50, 0, 5, "Radius (cm)"};
    AxisSpec cpaAxis = {50, 0.95, 1.0, "CPA"};
    AxisSpec tauAxis = {250, 0, 25, "Lifetime (cm)"};
    AxisSpec dcaAxis = {200, 0, 2, "DCA (cm)"};
    AxisSpec dcaxyAxis = {200, 0, 2, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {200, 0, 2, "DCA_{#it{z}} (cm)"};
    AxisSpec yAxis = {100, -1, 1, "Rapidity"};
    AxisSpec invMassAxisK0s = {800 / AxisConfig.cNbinsDiv, 0.46, 0.54, "Invariant Mass (GeV/#it{c}^2)"};    // K0s ~497.611
    AxisSpec invMassAxisChk892 = {900 / AxisConfig.cNbinsDiv, 0.5f, 1.4f, "Invariant Mass (GeV/#it{c}^2)"}; // chK(892) ~892
    AxisSpec invMassAxisReso = {1600 / AxisConfig.cNbinsDiv, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"};  // K1
    AxisSpec pidQAAxis = {130 / AxisConfig.cNbinsDivQA, -6.5, 6.5};
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

    // Secondary pion
    histos.add("QA/before/trkspionTPCPID", "TPC PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trkspionTOFPID", "TOF PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/before/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion 1 (positive) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});
    histos.add("QA/before/trkspionpT", "pT distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/before/trkspionDCAxy", "DCAxy distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/before/trkspionDCAz", "DCAz distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcazAxis});

    histos.add("QA/after/trkspionTPCPID", "TPC PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trkspionTOFPID", "TOF PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});
    histos.add("QA/after/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion 1 (positive) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});
    histos.add("QA/after/trkspionpT", "pT distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/after/trkspionDCAxy", "DCAxy distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcaxyAxis});
    histos.add("QA/after/trkspionDCAz", "DCAz distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcazAxis});

    // K0s pion 1
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

    // K0s pion 2
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
    histos.add("QA/before/hArmSecondary", "Armenteros distribution of secondary resonance", HistType::kTH2D, {AxisSpec{100, -1, 1, "alpha"}, {200, 0, 0.5, "qtArm"}});
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
    histos.add("QA/after/hArmSecondary", "Armenteros distribution of secondary resonance", HistType::kTH2D, {AxisSpec{100, -1, 1, "alpha"}, {200, 0, 0.5, "qtArm"}});
    histos.add("QA/after/hInvmassSecondary", "Invariant mass of unlike-sign secondary resonance", HistType::kTH1D, {invMassAxisK0s});

    // charged Kstar
    histos.add("QA/before/hpT_Kstar", "pT distribution of chK(892)", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/before/hy_Kstar", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
    histos.add("QA/before/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisChk892});

    histos.add("QA/after/hpT_Kstar", "pT distribution of chK(892)", HistType::kTH1D, {ptAxisQA});
    histos.add("QA/after/hy_Kstar", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
    histos.add("QA/after/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisChk892});

    // K1
    histos.add("QA/before/hpT_K1", "pT distribution of K1(1270)", HistType::kTH1F, {ptAxisQA});
    histos.add("QA/before/hy_K1", "Rapidity distribution of K1(1270)", HistType::kTH1F, {yAxis});
    histos.add("QA/before/K1CPA", "Cosine pointing angle distribution of K1(1270)", HistType::kTH1F, {cpaAxis});
    histos.add("QA/before/K1PtAsym", "pT asymmetry distribution of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QA/before/K1DalitzOS", "Dalitz plot of opposite-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/before/K1DalitzLS", "Dalitz plot of like-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/before/K1invmass", "Invariant mass of K1(1270) (US)", HistType::kTH1F, {invMassAxisReso});
    histos.add("QA/before/K1invmassLS", "Invariant mass of K1(1270) (LS)", HistType::kTH1F, {invMassAxisReso});

    histos.add("QA/after/hpT_K1", "pT distribution of K1(1270)", HistType::kTH1F, {ptAxisQA});
    histos.add("QA/after/hy_K1", "Rapidity distribution of K1(1270)", HistType::kTH1F, {yAxis});
    histos.add("QA/after/K1CPA", "Cosine pointing angle of K1(1270)", HistType::kTH1F, {cpaAxis});
    histos.add("QA/after/K1PtAsym", "pT asymmetry distribution of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QA/after/K1DalitzOS_Signal", "(Signal region) Dalitz plot of opposite-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/after/K1DalitzOS_SBLeft", "(SB left) Dalitz plot of opposite-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/after/K1DalitzOS_SBRight", "(SB right) Dalitz plot of opposite-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/after/K1DalitzLS_Signal", "(Signal region) Dalitz plot of like-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/after/K1DalitzLS_SBLeft", "(SB left) Dalitz plot of like-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/after/K1DalitzLS_SBRight", "(SB right) Dalitz plot of like-sign combination of  K1(1270)", HistType::kTH2F, {AxisSpec{300, 0.0, 3.0, "#it{M}^{2}_{K_{s}^{0}#pi_{sec}}"}, {300, 0.0, 3.0, "#it{M}^{2}_{#pi_{sec}#pi_{bach}}"}});
    histos.add("QA/after/K1invmass", "Invariant mass of K1(1270) (US)", HistType::kTH1F, {invMassAxisReso});
    histos.add("QA/after/K1invmassLS", "Invariant mass of K1(1270) (LS)", HistType::kTH1F, {invMassAxisReso});

    // Invariant mass
    histos.add("hInvmass_K1", "Invariant mass of K1(1270)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_Kstar", "Invariant mass of unlike-sign chK(892)", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisChk892});
    histos.add("hInvmass_K0s", "Invariant mass of unlike-sign K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s});

    ccdb->setURL(CCDBConfig.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in K1 Analysis Task";
    histos.print();
  } // init

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
    if (TrackCuts.cfgpTdepDCAzCut) {
      // Tuned on the LHC22f anchored MC LHC23d1d on primary pions. 7 Sigmas of the resolution
      if (std::abs(track.dcaZ()) > (0.004 + (0.013 / track.pt())))
        return false;
    } else {
      if (std::abs(track.dcaZ()) > TrackCuts.cfgMaxbDCAzToPVcut)
        return false;
    }
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
    //  if (candidate.pt() <= PIDCuts.cfgTOFMinPt)
    //    return true;

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

  // selection K0s
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

  K1MassRegion getK1MassRegion(float mass)
  {
    constexpr float SigMin = 1.22f;
    constexpr float SigMax = 1.32f;

    constexpr float SBLMin = 1.10f;
    constexpr float SBLMax = 1.18f;

    constexpr float SBRMin = 1.36f;
    constexpr float SBRMax = 1.44f;

    if (mass > SigMin && mass < SigMax)
      return K1MassRegion::Signal;
    if (mass > SBLMin && mass < SBLMax)
      return K1MassRegion::SBLeft;
    if (mass > SBRMin && mass < SBRMax)
      return K1MassRegion::SBRight;

    return K1MassRegion::Outside;
  }

  int count = 0;
  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType, typename TracksTypeK0s>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2, const TracksTypeK0s& dTracks3)
  {

    using TrackTarget = std::decay_t<TracksType>;

    histos.fill(HIST("QA/before/CentDist"), lCentrality);

    LorentzVectorSetXYZM lBachelor_pi, lDecayDaughter_K0s, lDecayDaughter_pi, lResoKstar, lResonanceK1, lDaughterRot, lResonanceRot;
    LorentzVectorSetXYZM lPairK0sPiSec, lPairK0sPiBach, lPairPiPi;

    std::vector<int> btrackIndices = {};
    std::vector<int> strackIndices = {};
    std::vector<int> k0sIndices = {};
    std::vector<int> chK892Indices = {};

    // check the existence of the pairs
    if (cfgHasPair && (dTracks1.size() < 1 || dTracks2.size() < 1 || dTracks3.size() < 1))
      return;

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
      btrackIndices.push_back(bTrack.index());
    } // bTrack

    for (const auto& sTrack : dTracks2) {
      auto trkspt = sTrack.pt();
      auto istrkshasTOF = sTrack.hasTOF();
      auto trksNSigmaPiTPC = sTrack.tpcNSigmaPi();
      auto trksNSigmaPiTOF = (istrkshasTOF) ? sTrack.tofNSigmaPi() : -999.;

      if constexpr (!IsMix) {
        if (cfgFillQAPlots) {
          // Secondary pion QA plots
          histos.fill(HIST("QA/before/trkspionTPCPID"), trkspt, trksNSigmaPiTPC);
          if (istrkshasTOF) {
            histos.fill(HIST("QA/before/trkspionTOFPID"), trkspt, trksNSigmaPiTOF);
            histos.fill(HIST("QA/before/trkspionTPCTOFPID"), trksNSigmaPiTPC, trksNSigmaPiTOF);
          }
          histos.fill(HIST("QA/before/trkspionpT"), trkspt);
          histos.fill(HIST("QA/before/trkspionDCAxy"), sTrack.dcaXY());
          histos.fill(HIST("QA/before/trkspionDCAz"), sTrack.dcaZ());
        }
      }

      if (!trackCut(sTrack))
        continue;
      if (!selectionPIDPion(sTrack))
        continue;

      if constexpr (!IsMix) {
        if (cfgFillQAPlots) {
          // Secondary pion QA plots after applying cuts
          histos.fill(HIST("QA/after/trkspionTPCPID"), trkspt, trksNSigmaPiTPC);
          if (istrkshasTOF) {
            histos.fill(HIST("QA/after/trkspionTOFPID"), trkspt, trksNSigmaPiTOF);
            histos.fill(HIST("QA/after/trkspionTPCTOFPID"), trksNSigmaPiTPC, trksNSigmaPiTOF);
          }
          histos.fill(HIST("QA/after/trkspionpT"), trkspt);
          histos.fill(HIST("QA/after/trkspionDCAxy"), sTrack.dcaXY());
          histos.fill(HIST("QA/after/trkspionDCAz"), sTrack.dcaZ());
        }
      }
      strackIndices.push_back(sTrack.index());
    } // sTrack

    for (const auto& k0sCand : dTracks3) {

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
          histos.fill(HIST("QA/before/hArmSecondary"), k0sCand.alpha(), k0sCand.qtarm());
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
          histos.fill(HIST("QA/after/hArmSecondary"), k0sCand.alpha(), k0sCand.qtarm());
          histos.fill(HIST("QA/after/hInvmassSecondary"), trkkMass);
        }
        histos.fill(HIST("hInvmass_K0s"), lCentrality, trkkpt, trkkMass);
      }
      k0sIndices.push_back(k0sCand.index());
    } // K0s

    for (const auto& btrackIndex : btrackIndices) {
      auto bTrack = dTracks1.rawIteratorAt(btrackIndex);
      for (const auto& strackIndex : strackIndices) {
        auto sTrack = dTracks2.rawIteratorAt(strackIndex);

        if (bTrack.index() == sTrack.index())
          continue;

        for (const auto& k0sIndex : k0sIndices) {
          auto k0sCand = dTracks3.rawIteratorAt(k0sIndex);

          auto posDauTrack = k0sCand.template posTrack_as<TrackTarget>();
          auto negDauTrack = k0sCand.template negTrack_as<TrackTarget>();

          if (bTrack.index() == posDauTrack.index() || bTrack.index() == negDauTrack.index())
            continue;
          if (sTrack.index() == posDauTrack.index() || sTrack.index() == negDauTrack.index())
            continue;

          lBachelor_pi = LorentzVectorSetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), MassPionCharged);
          lDecayDaughter_pi = LorentzVectorSetXYZM(sTrack.px(), sTrack.py(), sTrack.pz(), MassPionCharged);
          lDecayDaughter_K0s = LorentzVectorSetXYZM(k0sCand.px(), k0sCand.py(), k0sCand.pz(), MassK0Short);
          lResoKstar = lDecayDaughter_pi + lDecayDaughter_K0s;
          lResonanceK1 = lResoKstar + lBachelor_pi;

          // QA plots for Kstar
          if constexpr (!IsMix) {
            if (cfgFillQAPlots) {
              histos.fill(HIST("QA/before/hpT_Kstar"), lResoKstar.Pt());
              histos.fill(HIST("QA/before/hy_Kstar"), lResoKstar.Rapidity());
              histos.fill(HIST("QA/before/kstarinvmass"), lResoKstar.M());
            }
          }

          if (lResoKstar.Rapidity() > KstarCuts.cfgKstarMaxRap || lResoKstar.Rapidity() < KstarCuts.cfgKstarMinRap)
            continue;
          if (lResoKstar.Pt() < KstarCuts.cfgKstarMinPt)
            continue;
          if (std::fabs(lResoKstar.M() - MassKPlusStar892) > KstarCuts.cfgKstarMassWindow)
            continue;

          if constexpr (!IsMix) {
            if (cfgFillQAPlots) {
              histos.fill(HIST("QA/after/hpT_Kstar"), lResoKstar.Pt());
              histos.fill(HIST("QA/after/hy_Kstar"), lResoKstar.Rapidity());
              histos.fill(HIST("QA/after/kstarinvmass"), lResoKstar.M());
            }
            histos.fill(HIST("hInvmass_Kstar"), lCentrality, lResoKstar.Pt(), lResoKstar.M());
          } // IsMix

          lPairK0sPiSec = lDecayDaughter_K0s + lDecayDaughter_pi;
          lPairK0sPiBach = lDecayDaughter_K0s + lBachelor_pi;
          lPairPiPi = lDecayDaughter_pi + lBachelor_pi;

          float m2K0sPiSec = lPairK0sPiSec.M2();
          // float m2K0sPiBach = lPairK0sPiBach.M2();
          float m2PiPi = lPairPiPi.M2();

          // QA plots for K1
          if constexpr (!IsMix) {
            if (cfgFillQAPlots) {
              histos.fill(HIST("QA/before/hpT_K1"), lResonanceK1.Pt());
              histos.fill(HIST("QA/before/hy_K1"), lResonanceK1.Rapidity());
              if (bTrack.sign() * sTrack.sign() < 0) {
                histos.fill(HIST("QA/before/K1invmass"), lResonanceK1.M());
                histos.fill(HIST("QA/before/K1DalitzOS"), m2K0sPiSec, m2PiPi);
              } else {
                histos.fill(HIST("QA/before/K1invmassLS"), lResonanceK1.M());
                histos.fill(HIST("QA/before/K1DalitzLS"), m2K0sPiSec, m2PiPi);
              }
            }
          } // IsMix

          if (lResonanceK1.Rapidity() > K1Cuts.cfgK1MaxRap || lResonanceK1.Rapidity() < K1Cuts.cfgK1MinRap)
            continue;
          if (lResonanceK1.Pt() < K1Cuts.cfgK1MinPt)
            continue;
          if (lPairPiPi.Pt() < cfgPiPiMinPt)
            continue;

          auto k1Region = getK1MassRegion(lResonanceK1.M());

          if constexpr (!IsMix) {
            if (cfgFillQAPlots) {
              histos.fill(HIST("QA/after/hpT_K1"), lResonanceK1.Pt());
              histos.fill(HIST("QA/after/hy_K1"), lResonanceK1.Rapidity());
              if (bTrack.sign() * sTrack.sign() < 0) {
                histos.fill(HIST("QA/after/K1invmass"), lResonanceK1.M());
                // histos.fill(HIST("QA/after/K1DalitzOS"), m2K0sPiSec, m2PiPi);
                if (k1Region == K1MassRegion::Signal) {
                  histos.fill(HIST("QA/after/K1DalitzOS_Signal"), m2K0sPiSec, m2PiPi);
                } else if (k1Region == K1MassRegion::SBLeft) {
                  histos.fill(HIST("QA/after/K1DalitzOS_SBLeft"), m2K0sPiSec, m2PiPi);
                } else if (k1Region == K1MassRegion::SBRight) {
                  histos.fill(HIST("QA/after/K1DalitzOS_SBRight"), m2K0sPiSec, m2PiPi);
                }
              } else {
                histos.fill(HIST("QA/after/K1invmassLS"), lResonanceK1.M());
                // histos.fill(HIST("QA/after/K1DalitzLS"), m2K0sPiSec, m2PiPi);
                if (k1Region == K1MassRegion::Signal) {
                  histos.fill(HIST("QA/after/K1DalitzLS_Signal"), m2K0sPiSec, m2PiPi);
                } else if (k1Region == K1MassRegion::SBLeft) {
                  histos.fill(HIST("QA/after/K1DalitzLS_SBLeft"), m2K0sPiSec, m2PiPi);
                } else if (k1Region == K1MassRegion::SBRight) {
                  histos.fill(HIST("QA/after/K1DalitzLS_SBRight"), m2K0sPiSec, m2PiPi);
                }
              }
            }

            if (bTrack.sign() * sTrack.sign() < 0) {
              // bTrack sign minus for particle, plus for anti-particle
              unsigned int typeK1 = bTrack.sign() < 0 ? BinType::kK1P : BinType::kK1A;
              histos.fill(HIST("hInvmass_K1"), typeK1, lCentrality, lResonanceK1.Pt(), lResonanceK1.M());
            } else {
              unsigned int typeK1 = bTrack.sign() < 0 ? BinType::kK1P_Like : BinType::kK1A_Like;
              histos.fill(HIST("hInvmass_K1"), typeK1, lCentrality, lResonanceK1.Pt(), lResonanceK1.M());
            }

            if (BkgEstimationConfig.cfgFillRotBkg) {
              for (int i = 0; i < BkgEstimationConfig.cfgNrotBkg; i++) {
                auto lRotAngle = BkgEstimationConfig.cfgMinRot + i * ((BkgEstimationConfig.cfgMaxRot - BkgEstimationConfig.cfgMinRot) / (BkgEstimationConfig.cfgNrotBkg - 1));
                if (cfgFillQAPlots) {
                  histos.fill(HIST("QA/RotBkg/hRotBkg"), lRotAngle);
                }
                if (BkgEstimationConfig.cfgRotPion) {
                  lDaughterRot = lBachelor_pi;
                  ROOT::Math::RotationZ rot(lRotAngle);
                  auto p3 = rot * lDaughterRot.Vect();
                  lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                  lResonanceRot = lDaughterRot + lResoKstar;
                } else {
                  lDaughterRot = lResoKstar;
                  ROOT::Math::RotationZ rot(lRotAngle);
                  auto p3 = rot * lDaughterRot;
                  lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                  lResonanceRot = lBachelor_pi + lDaughterRot;
                }
                unsigned int typeK1 = bTrack.sign() < 0 ? BinType::kK1P_Rot : BinType::kK1A_Rot;
                histos.fill(HIST("hInvmass_K1"), typeK1, lCentrality, lResonanceRot.Pt(), lResonanceRot.M());

              } // NrotBkg
            } // cfgFillRotBkg
          } // IsMix
        } // k0sIndex
      } // strackIndex
    } // btrackIndex

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

    fillHistograms<false, false>(collision, tracks, tracks, v0s);
  }
  PROCESS_SWITCH(K1analysis, processData, "Process Event for data without Partitioning", true);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K1analysis>(cfgc)};
}
