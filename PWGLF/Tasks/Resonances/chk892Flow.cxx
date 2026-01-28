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
/// \file chk892Flow.cxx
/// \brief Reconstruction of track-track decay resonance candidates
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <Bong-Hwi.Lim@cern.ch>
///

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/collisionCuts.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
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
#include "Framework/EndOfStreamContext.h"
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
#include "TRandom3.h"
#include "TVector2.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>

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

struct Chk892Flow {
  enum BinType : unsigned int {
    kKstarP = 0,
    kKstarN,
    kKstarP_Mix,
    kKstarN_Mix,
    kKstarP_Rot,
    kKstarN_Rot,
    kTYEnd
  };

  enum class K0sCut {
    DauDCA,        // lDauDCA <= cfgSecondaryDauDCAMax
    PosDCAtoPVMin, // lDauPosDCAtoPV >= min
    NegDCAtoPVMin, // lDauNegDCAtoPV >= min
    RadiusWindow,  // cfgSecondaryRadiusMin <= lRadius <= cfgSecondaryRadiusMax
    DCAtoPVMax,    // lDCAtoPV <= max
    CPAMin,        // lCPA >= min
    ProperTauMax,  // lPropTauK0s <= max
    Armenteros,    // qtarm >= param * |alpha|
    MassWindow,    // |mK0s - m0| <= window
    LambdaMassHypo // NOT(lambda-window)
  };

  std::array<std::shared_ptr<TH1>, 10> hN1NoCut{};
  std::array<std::shared_ptr<TH1>, 10> hN1Pass{};

  static constexpr const char* cutTag[] = {
    "DauDCA", "PosDCA", "NegDCA", "Radius", "DCAtoPV", "CPA", "Tau", "Arm", "Mass", "LambdaHypo"};
  static constexpr K0sCut kCutsToTest[] = {
    K0sCut::DauDCA, K0sCut::PosDCAtoPVMin, K0sCut::NegDCAtoPVMin,
    K0sCut::RadiusWindow, K0sCut::DCAtoPVMax, K0sCut::CPAMin,
    K0sCut::ProperTauMax, K0sCut::Armenteros, K0sCut::MassWindow, K0sCut::LambdaMassHypo};

  static constexpr size_t NCuts = sizeof(kCutsToTest) / sizeof(kCutsToTest[0]);

  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors>;
  // using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi>;
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
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::framework::O2DatabasePDG> pdg;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  } CCDBConfig;
  // Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  // Configurables
  struct : ConfigurableGroup {
    ConfigurableAxis cfgBinsPt{"cfgBinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsPtQA{"cfgBinsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsCent{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
    ConfigurableAxis cfgBinsVtxZ{"cfgBinsVtxZ", {VARIABLE_WIDTH, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "Binning of the z-vertex axis"};
    ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0.0, 3.00065, 4.28798, 6.14552, 7.6196, 8.90942, 10.0897, 11.2002, 12.2709, 13.3167, 14.4173, 23.2518}, "Binning of the impact parameter axis"};
    ConfigurableAxis cfgBinsOccu{"cfgBinsOccu", {VARIABLE_WIDTH, 0, 500, 1000, 2500, 10000, 999999}, "Binning of the occupancy axis"};
    Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};
    Configurable<int> cNbinsDivQA{"cNbinsDivQA", 1, "Integer to divide the number of bins for QA"};
    ConfigurableAxis cfgAxisV2{"cfgAxisV2", {200, -1, 1}, "Binning of the v2 axis (+-1 for EP method)"};
    ConfigurableAxis cfgAxisPhi{"cfgAxisPhi", {8, 0, constants::math::PI}, "Binning of the #phi axis"};
  } AxisConfig;

  struct : ConfigurableGroup {
    Configurable<bool> cfgFillQAPlots{"cfgFillQAPlots", true, "Fill QA plots"};
    Configurable<bool> cfgQvecSel{"cfgQvecSel", true, "Reject events when no QVector"};
    Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};
    Configurable<bool> cfgFillAdditionalAxis{"cfgFillAdditionalAxis", false, "Fill additional axis"};
    Configurable<bool> cfgUseScalProduct{"cfgUseScalProduct", false, "Use scalar product method"};
    Configurable<bool> cfgFillOccupancy{"cfgFillOccupancy", false, "Fill Occupancy"};
  } AnalysisConfig;

  // Event cuts
  o2::analysis::CollisonCuts colCuts;
  struct : ConfigurableGroup {
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<float> cfgEvtZvtxMC{"cfgEvtZvtxMC", 10.f, "MC Evt sel: Max z-vertex (cm)"};
    Configurable<bool> cfgIsPhysicalPrimaryMC{"cfgIsPhysicalPrimaryMC", true, "Physical primary selection for MC parents"};
    //    Configurable<int> cfgEvtOccupancyInTimeRangeMax{"cfgEvtOccupancyInTimeRangeMax", -1, "Evt sel: maximum track occupancy"};
    //    Configurable<int> cfgEvtOccupancyInTimeRangeMin{"cfgEvtOccupancyInTimeRangeMin", -1, "Evt sel: minimum track occupancy"};
    Configurable<bool> cfgEvtTriggerCheck{"cfgEvtTriggerCheck", false, "Evt sel: check for trigger"};
    Configurable<bool> cfgEvtOfflineCheck{"cfgEvtOfflineCheck", true, "Evt sel: check for offline selection"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", false, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtTFBorderCut{"cfgEvtTFBorderCut", false, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtUseITSTPCvertex{"cfgEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
    Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", true, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", true, "Evt sel: apply pileup rejection"};
    Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtCollInTimeRangeStandard{"cfgEvtCollInTimeRangeStandard", true, "Evt sel: apply NoCollInTimeRangeStandard"};
    Configurable<float> cfgEventCentralityMin{"cfgEventCentralityMin", 0.0f, "Event sel: minimum centrality"};
    Configurable<float> cfgEventCentralityMax{"cfgEventCentralityMax", 80.0f, "Event sel: maximum centrality"};
    Configurable<bool> cfgEvtUseRCTFlagChecker{"cfgEvtUseRCTFlagChecker", false, "Evt sel: use RCT flag checker"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } EventCuts;
  RCTFlagsChecker rctChecker;

  /// PID Selections, pion
  struct : ConfigurableGroup {
    Configurable<bool> cfgTPConly{"cfgTPConly", false, "Use only TPC for PID"};                                     // bool
    Configurable<float> cfgMaxTPCnSigmaPion{"cfgMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};                 // TPC
    Configurable<float> cfgMaxTOFnSigmaPion{"cfgMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                 // TOF
    Configurable<float> cfgNsigmaCutCombinedPion{"cfgNsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"}; // Combined
    Configurable<bool> cfgTOFVeto{"cfgTOFVeto", true, "TOF Veto, if false, TOF is nessessary for PID selection"};   // TOF Veto
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
    Configurable<bool> cfgpTdepDCAxyCut{"cfgpTdepDCAxyCut", false, "pT-dependent DCAxy cut"};
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
    // DCA to PV
    Configurable<float> cfgMaxbDCArToPVcut{"cfgMaxbDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
    Configurable<float> cfgMaxbDCAzToPVcut{"cfgMaxbDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
  } TrackCuts;

  // Secondary Selection
  struct : ConfigurableGroup {
    Configurable<bool> cfgReturnFlag{"cfgReturnFlag", false, "Return Flag for debugging"};
    Configurable<bool> cfgSecondaryRequire{"cfgSecondaryRequire", true, "Secondary cuts on/off"};
    Configurable<bool> cfgSecondaryArmenterosCut{"cfgSecondaryArmenterosCut", true, "cut on Armenteros-Podolanski graph"};
    Configurable<bool> cfgSecondaryCrossMassHypothesisCut{"cfgSecondaryCrossMassHypothesisCut", false, "Apply cut based on the lambda mass hypothesis"};

    Configurable<bool> cfgByPassDauPIDSelection{"cfgByPassDauPIDSelection", true, "Bypass Daughters PID selection"};
    Configurable<bool> cfgByPassDauRapiditySelection{"cfgByPassDauRapiditySelection", false, "Bypass Daughters Rapidity selection"};
    Configurable<float> cfgSecondaryDauDCAMax{"cfgSecondaryDauDCAMax", 0.2, "Maximum DCA Secondary daughters to PV"};
    Configurable<float> cfgSecondaryDauPosDCAtoPVMin{"cfgSecondaryDauPosDCAtoPVMin", 0.1, "Minimum DCA Secondary positive daughters to PV"};
    Configurable<float> cfgSecondaryDauNegDCAtoPVMin{"cfgSecondaryDauNegDCAtoPVMin", 0.1, "Minimum DCA Secondary negative daughters to PV"};

    Configurable<float> cfgSecondaryPtMin{"cfgSecondaryPtMin", 0.f, "Minimum transverse momentum of Secondary"};
    Configurable<float> cfgSecondaryRapidityMax{"cfgSecondaryRapidityMax", 0.8, "Maximum rapidity of Secondary"};
    Configurable<float> cfgSecondaryDauRapidityMax{"cfgSecondaryDauRapidityMax", 0.3, "Maximum rapidity of Secondary daughters"};
    Configurable<float> cfgSecondaryRadiusMin{"cfgSecondaryRadiusMin", 0.0, "Minimum transverse radius of Secondary"};
    Configurable<float> cfgSecondaryRadiusMax{"cfgSecondaryRadiusMax", 999.9, "Maximum transverse radius of Secondary"};
    Configurable<float> cfgSecondaryCosPAMin{"cfgSecondaryCosPAMin", 0.995, "Mininum cosine pointing angle of Secondary"};
    Configurable<float> cfgSecondaryDCAtoPVMax{"cfgSecondaryDCAtoPVMax", 0.4, "Maximum DCA Secondary to PV"};
    Configurable<float> cfgSecondaryProperLifetimeMax{"cfgSecondaryProperLifetimeMax", 20., "Maximum Secondary Lifetime"};
    Configurable<float> cfgSecondaryparamArmenterosCut{"cfgSecondaryparamArmenterosCut", 0.2, "parameter for Armenteros Cut"};
    Configurable<float> cfgSecondaryMassWindow{"cfgSecondaryMassWindow", 0.03, "Secondary inv mass selection window"};
    Configurable<float> cfgSecondaryCrossMassCutWindow{"cfgSecondaryCrossMassCutWindow", 0.05, "Secondary inv mass selection window with (anti)lambda hypothesis"};
  } SecondaryCuts;

  // K* selection
  struct : ConfigurableGroup {
    Configurable<float> cfgKstarMaxRap{"cfgKstarMaxRap", 0.5, "Kstar maximum rapidity"};
    Configurable<float> cfgKstarMinRap{"cfgKstarMinRap", -0.5, "Kstar minimum rapidity"};
  } KstarCuts;

  // Confs from flow analysis
  struct : ConfigurableGroup {
    Configurable<int> cfgnMods{"cfgnMods", 2, "The number of modulations of interest starting from 2"};
    Configurable<int> cfgNQvec{"cfgNQvec", 7, "The number of total Qvectors for looping over the task"};

    Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0C", "The name of detector to be analyzed"};
    Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "The name of detector for reference A"};
    Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "The name of detector for reference B"};
  } EventPlaneConfig;

  // Bkg estimation
  struct : ConfigurableGroup {
    Configurable<bool> cfgFillRotBkg{"cfgFillRotBkg", true, "Fill rotated background"};
    Configurable<float> cfgMinRot{"cfgMinRot", 5.0 * constants::math::PI / 6.0, "Minimum of rotation"};
    Configurable<float> cfgMaxRot{"cfgMaxRot", 7.0 * constants::math::PI / 6.0, "Maximum of rotation"};
    Configurable<bool> cfgRotPion{"cfgRotPion", true, "Rotate pion"};
    Configurable<int> cfgNrotBkg{"cfgNrotBkg", 9, "Number of rotated copies (background) per each original candidate"};
  } BkgEstimationConfig;

  int lDetId;
  int lRefAId;
  int lRefBId;

  int lQvecDetInd;
  int lQvecRefAInd;
  int lQvecRefBInd;

  float lCentrality;

  // PDG code
  int kPDGK0s = kK0Short;
  int kPDGK0 = kK0;
  int kKstarPlus = o2::constants::physics::Pdg::kKPlusStar892;

  void init(o2::framework::InitContext&)
  {
    lCentrality = -999;

    colCuts.setCuts(EventCuts.cfgEvtZvtx, EventCuts.cfgEvtTriggerCheck, EventCuts.cfgEvtOfflineCheck, /*checkRun3*/ true, /*triggerTVXsel*/ false, /*EventCuts.cfgEvtOccupancyInTimeRangeMax*/ false, /*EventCuts.cfgEvtOccupancyInTimeRangeMin*/ false);
    colCuts.init(&histos);
    colCuts.setTriggerTVX(EventCuts.cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(EventCuts.cfgEvtTFBorderCut);
    colCuts.setApplyITSTPCvertex(EventCuts.cfgEvtUseITSTPCvertex);
    colCuts.setApplyZvertexTimedifference(EventCuts.cfgEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(EventCuts.cfgEvtPileupRejection);
    colCuts.setApplyNoITSROBorderCut(EventCuts.cfgEvtNoITSROBorderCut);
    colCuts.setApplyCollInTimeRangeStandard(EventCuts.cfgEvtCollInTimeRangeStandard);
    colCuts.printCuts();

    rctChecker.init(EventCuts.cfgEvtRCTFlagCheckerLabel, EventCuts.cfgEvtRCTFlagCheckerZDCCheck, EventCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    AxisSpec centAxis = {AxisConfig.cfgBinsCent, "T0M (%)"};
    AxisSpec vtxzAxis = {AxisConfig.cfgBinsVtxZ, "Z Vertex (cm)"};
    AxisSpec epAxis = {100, -1.0 * constants::math::PI, constants::math::PI};
    AxisSpec ptAxis = {AxisConfig.cfgBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {AxisConfig.cfgBinsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec v2Axis = {AxisConfig.cfgAxisV2, "#v_{2}"};
    AxisSpec phiAxis = {AxisConfig.cfgAxisPhi, "2(#phi-#Psi_{2})"};
    AxisSpec occuAxis = {AxisConfig.cfgBinsOccu, "Occupancy"};
    AxisSpec radiusAxis = {50, 0, 5, "Radius (cm)"};
    AxisSpec cpaAxis = {30, 0.97, 1.0, "CPA"};
    AxisSpec tauAxis = {250, 0, 25, "Lifetime (cm)"};
    AxisSpec dcaAxis = {100, 0, 2, "DCA (cm)"};
    AxisSpec dcaxyAxis = {100, 0, 1, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {200, 0, 2, "DCA_{#it{z}} (cm)"};
    AxisSpec yAxis = {50, -1, 1, "Rapidity"};
    AxisSpec invMassAxisK0s = {800 / AxisConfig.cNbinsDiv, 0.46, 0.54, "Invariant Mass (GeV/#it{c}^2)"};  // K0s ~497.611
    AxisSpec invMassAxisReso = {900 / AxisConfig.cNbinsDiv, 0.5f, 1.4f, "Invariant Mass (GeV/#it{c}^2)"}; // chK(892) ~892
    AxisSpec pidQAAxis = {130 / AxisConfig.cNbinsDivQA, -6.5, 6.5};

    // THnSparse
    AxisSpec axisType = {BinType::kTYEnd, 0, BinType::kTYEnd, "Type of bin with charge and mix"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};

    if (SecondaryCuts.cfgReturnFlag) {
      histos.add("QA/K0sCutCheck", "Check K0s cut", HistType::kTH1D, {AxisSpec{13, -0.5, 12.5, "Check"}});
    }
    histos.add("QA/before/CentDist", "Centrality distribution", {HistType::kTH1D, {centAxis}});
    histos.add("QA/before/VtxZ", "z-vertex distribution", {HistType::kTH1D, {vtxzAxis}});
    histos.add("QA/before/Occupancy", "Occupancy distribution", {HistType::kTH1D, {occuAxis}});

    // EventPlane
    histos.add("QA/EP/hEPDet", "Event plane distribution of FT0C (Det = A)", {HistType::kTH2D, {centAxis, epAxis}});
    histos.add("QA/EP/hEPB", "Event plane distribution of TPCpos (B)", {HistType::kTH2D, {centAxis, epAxis}});
    histos.add("QA/EP/hEPC", "Event plane distribution of TPCneg (C)", {HistType::kTH2D, {centAxis, epAxis}});
    histos.add("QA/EP/hEPResAB", "cos(n(A-B))", {HistType::kTH2D, {centAxis, epAxis}});
    histos.add("QA/EP/hEPResAC", "cos(n(A-C))", {HistType::kTH2D, {centAxis, epAxis}});
    histos.add("QA/EP/hEPResBC", "cos(n(B-C))", {HistType::kTH2D, {centAxis, epAxis}});

    if (AnalysisConfig.cfgUseScalProduct) {
      histos.add("QA/EP/hEPSPResAB", "cos(n(A-B))", {HistType::kTH2D, {centAxis, epAxis}});
      histos.add("QA/EP/hEPSPResAC", "cos(n(A-C))", {HistType::kTH2D, {centAxis, epAxis}});
      histos.add("QA/EP/hEPSPResBC", "cos(n(B-C))", {HistType::kTH2D, {centAxis, epAxis}});
    }

    if (AnalysisConfig.cfgFillQAPlots) {
      // Rotated background
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
      histos.add("QA/before/hy_Secondary", "Rapidity distribution of secondary resonance", HistType::kTH1D, {yAxis});
      histos.add("QA/before/hCPASecondary", "Cosine pointing angle distribution of secondary resonance", HistType::kTH1D, {cpaAxis});
      histos.add("QA/before/hDCAtoPVSecondary", "DCA to PV distribution of secondary resonance", HistType::kTH1D, {dcaAxis});
      histos.add("QA/before/hPropTauSecondary", "Proper Lifetime distribution of secondary resonance", HistType::kTH1D, {tauAxis});
      histos.add("QA/before/hInvmassSecondary", "Invariant mass of unlike-sign secondary resonance", HistType::kTH1D, {invMassAxisK0s});

      histos.add("QA/after/hDauDCASecondary", "DCA of daughters of secondary resonance", HistType::kTH1D, {dcaAxis});
      histos.add("QA/after/hy_Secondary", "Rapidity distribution of secondary resonance", HistType::kTH1D, {yAxis});
      histos.add("QA/after/hCPASecondary", "Cosine pointing angle distribution of secondary resonance", HistType::kTH1D, {cpaAxis});
      histos.add("QA/after/hDCAtoPVSecondary", "DCA to PV distribution of secondary resonance", HistType::kTH1D, {dcaAxis});
      histos.add("QA/after/hPropTauSecondary", "Proper Lifetime distribution of secondary resonance", HistType::kTH1D, {tauAxis});
      histos.add("QA/after/hInvmassSecondary", "Invariant mass of unlike-sign secondary resonance", HistType::kTH1D, {invMassAxisK0s});

      // Mass QA (quick check)
      histos.add("QA/before/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
      histos.add("QA/before/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
      histos.add("QA/before/k0sv2vsinvmass", "Invariant mass vs v2 of unlike-sign K0s", HistType::kTH2D, {invMassAxisK0s, v2Axis});
      histos.add("QA/before/kstarv2vsinvmass", "Invariant mass vs v2 of unlike-sign chK(892)", HistType::kTH2D, {invMassAxisReso, v2Axis});
      histos.add("QA/before/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

      histos.add("QA/after/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
      histos.add("QA/after/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
      histos.add("QA/after/k0sv2vsinvmass", "Invariant mass vs v2 of unlike-sign K0s", HistType::kTH2D, {invMassAxisK0s, v2Axis});
      histos.add("QA/after/kstarv2vsinvmass", "Invariant mass vs v2 of unlike-sign chK(892)", HistType::kTH2D, {invMassAxisReso, v2Axis});
      histos.add("QA/after/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

      LOG(info) << "Size of the histograms in spectraTOF";
      histos.print();
    }

    // Invariant mass nSparse
    if (AnalysisConfig.cfgFillAdditionalAxis) {
      if (AnalysisConfig.cfgFillOccupancy) {
        histos.add("hInvmass_Kstar", "Invariant mass of unlike-sign chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis, phiAxis, occuAxis});
        histos.add("hInvmass_K0s", "Invariant mass of unlike-sign K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s, v2Axis, phiAxis, occuAxis});
        if (doprocessMC) {
          histos.add("hInvmass_Kstar_MC", "Invariant mass of unlike chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis, phiAxis, occuAxis});
          histos.add("hInvmass_K0s_MC", "Invariant mass of unlike K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso, v2Axis, phiAxis, occuAxis});
        }
      } else {
        histos.add("hInvmass_Kstar", "Invariant mass of unlike-sign chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis, phiAxis});
        histos.add("hInvmass_K0s", "Invariant mass of unlike-sign K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s, v2Axis, phiAxis});
        if (doprocessMC) {
          histos.add("hInvmass_Kstar_MC", "Invariant mass of unlike chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis, phiAxis});
          histos.add("hInvmass_K0s_MC", "Invariant mass of unlike K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso, v2Axis, phiAxis});
        }
      }
    } else {
      if (AnalysisConfig.cfgFillOccupancy) {
        histos.add("hInvmass_Kstar", "Invariant mass of unlike-sign chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis, occuAxis});
        histos.add("hInvmass_K0s", "Invariant mass of unlike-sign K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s, v2Axis, occuAxis});
        if (doprocessMC) {
          histos.add("hInvmass_Kstar_MC", "Invariant mass of unlike chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis, occuAxis});
          histos.add("hInvmass_K0s_MC", "Invariant mass of unlike K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s, v2Axis, occuAxis});
        }
      } else {
        histos.add("hInvmass_Kstar", "Invariant mass of unlike-sign chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis});
        histos.add("hInvmass_K0s", "Invariant mass of unlike-sign K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s, v2Axis});
        if (doprocessMC) {
          histos.add("hInvmass_Kstar_MC", "Invariant mass of unlike chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso, v2Axis});
          histos.add("hInvmass_K0s_MC", "Invariant mass of unlike K0s", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisK0s, v2Axis});
        }
      }
    }

    lDetId = getlDetId(EventPlaneConfig.cfgQvecDetName);
    lRefAId = getlDetId(EventPlaneConfig.cfgQvecRefAName);
    lRefBId = getlDetId(EventPlaneConfig.cfgQvecRefBName);

    if (lDetId == lRefAId || lDetId == lRefBId || lRefAId == lRefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      lDetId = 0;
      lRefAId = 4;
      lRefBId = 5;
    }
    if (EventPlaneConfig.cfgNQvec < EventPlaneConfig.cfgnMods) {
      LOG(fatal) << "nMode must be larger than 1, current input (cfgNQvec): " << EventPlaneConfig.cfgNQvec;
    }
    LOGF(info, "lDetId: %d, lRefAId: %d, lRefBId: %d", lDetId, lRefAId, lRefBId);

    ccdb->setURL(CCDBConfig.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in chK(892) Analysis Task";
    histos.print();
  } // init

  const int kCentFT0C = 1;
  const int kCentFT0M = 2;
  const float kInvalidCentrality = -999.f;

  template <typename CollisionType>
  float getCentrality(CollisionType const& collision)
  {
    if (AnalysisConfig.cfgCentEst == kCentFT0C) {
      return collision.centFT0C();
    } else if (AnalysisConfig.cfgCentEst == kCentFT0M) {
      return collision.centFT0M();
    } else {
      return kInvalidCentrality;
    }
  }

  template <typename DetNameType>
  int getlDetId(DetNameType const& name)
  {
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos") {
      return 4;
    } else if (name.value == "TPCneg") {
      return 5;
    } else {
      return false;
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
    auto lDauDCA = candidate.dcaV0daughters();
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
      if (lDauDCA > SecondaryCuts.cfgSecondaryDauDCAMax)
        return false;
      if (lDauPosDCAtoPV < SecondaryCuts.cfgSecondaryDauPosDCAtoPVMin)
        return false;
      if (lDauNegDCAtoPV < SecondaryCuts.cfgSecondaryDauNegDCAtoPVMin)
        return false;
      if (lPt < SecondaryCuts.cfgSecondaryPtMin)
        return false;
      if (std::fabs(lRapidity) > SecondaryCuts.cfgSecondaryRapidityMax)
        return false;
      if (lRadius < SecondaryCuts.cfgSecondaryRadiusMin || lRadius > SecondaryCuts.cfgSecondaryRadiusMax)
        return false;
      if (lDCAtoPV > SecondaryCuts.cfgSecondaryDCAtoPVMax)
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
      if (candidate.qtarm() < SecondaryCuts.cfgSecondaryparamArmenterosCut * std::abs(candidate.alpha())) {
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
  bool isTrueKstar(const TrackTemplate& bTrack, const V0Template& k0sCand)
  {
    if (std::abs(bTrack.pdgCode()) != kPiPlus) // Are you pion?
      return false;
    if (std::abs(k0sCand.pdgCode()) != kPDGK0s) // Are you K0s?
      return false;

    auto motherbTrack = bTrack.template mothers_as<aod::McParticles>();
    auto motherkV0 = k0sCand.template mothers_as<aod::McParticles>();

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

  int count = 0;

  template <bool IsMC, bool IsMCTrue, bool IsMix, typename CollisionType, typename TracksType, typename TracksTypeK0s>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksTypeK0s& dTracks2, int nmode)
  {
    histos.fill(HIST("QA/before/CentDist"), lCentrality);
    histos.fill(HIST("QA/before/Occupancy"), collision.trackOccupancyInTimeRange());

    lQvecDetInd = lDetId * 4 + 3 + (nmode - 2) * EventPlaneConfig.cfgNQvec * 4;
    lQvecRefAInd = lRefAId * 4 + 3 + (nmode - 2) * EventPlaneConfig.cfgNQvec * 4;
    lQvecRefBInd = lRefBId * 4 + 3 + (nmode - 2) * EventPlaneConfig.cfgNQvec * 4;

    double lEPDet = std::atan2(collision.qvecIm()[lQvecDetInd], collision.qvecRe()[lQvecDetInd]) / static_cast<float>(nmode);
    double lEPRefB = std::atan2(collision.qvecIm()[lQvecRefAInd], collision.qvecRe()[lQvecRefAInd]) / static_cast<float>(nmode);
    double lEPRefC = std::atan2(collision.qvecIm()[lQvecRefBInd], collision.qvecRe()[lQvecRefBInd]) / static_cast<float>(nmode);

    double lEPResAB = std::cos(static_cast<float>(nmode) * (lEPDet - lEPRefB));
    double lEPResAC = std::cos(static_cast<float>(nmode) * (lEPDet - lEPRefC));
    double lEPResBC = std::cos(static_cast<float>(nmode) * (lEPRefB - lEPRefC));

    // EP method
    histos.fill(HIST("QA/EP/hEPDet"), lCentrality, lEPDet);
    histos.fill(HIST("QA/EP/hEPB"), lCentrality, lEPRefB);
    histos.fill(HIST("QA/EP/hEPC"), lCentrality, lEPRefC);
    histos.fill(HIST("QA/EP/hEPResAB"), lCentrality, lEPResAB);
    histos.fill(HIST("QA/EP/hEPResAC"), lCentrality, lEPResAC);
    histos.fill(HIST("QA/EP/hEPResBC"), lCentrality, lEPResBC);
    // Scalar product method
    if (AnalysisConfig.cfgUseScalProduct) {
      double lEPSPResAB = (collision.qvecRe()[lQvecDetInd] * collision.qvecRe()[lQvecRefAInd] + collision.qvecIm()[lQvecDetInd] * collision.qvecIm()[lQvecRefAInd]);
      double lEPSPResAC = (collision.qvecRe()[lQvecDetInd] * collision.qvecRe()[lQvecRefBInd] + collision.qvecIm()[lQvecDetInd] * collision.qvecIm()[lQvecRefBInd]);
      double lEPSPResBC = (collision.qvecRe()[lQvecRefAInd] * collision.qvecRe()[lQvecRefBInd] + collision.qvecIm()[lQvecRefAInd] * collision.qvecIm()[lQvecRefBInd]);

      histos.fill(HIST("QA/EP/hEPSPResAB"), lCentrality, lEPSPResAB);
      histos.fill(HIST("QA/EP/hEPSPResAC"), lCentrality, lEPSPResAC);
      histos.fill(HIST("QA/EP/hEPSPResBC"), lCentrality, lEPSPResBC);
    }

    LorentzVectorSetXYZM lDecayDaughter1, lDecayDaughter2, lResoSecondary, lDecayDaughter_bach, lResoKstar, lDaughterRot, lResonanceRot;
    std::vector<int> trackIndicies = {};
    std::vector<int> k0sIndicies = {};

    for (const auto& bTrack : dTracks1) {
      auto trkbpt = bTrack.pt();
      auto istrkbhasTOF = bTrack.hasTOF();
      auto trkbNSigmaPiTPC = bTrack.tpcNSigmaPi();
      auto trkbNSigmaPiTOF = (istrkbhasTOF) ? bTrack.tofNSigmaPi() : -999.;

      if constexpr (!IsMix) {
        if (AnalysisConfig.cfgFillQAPlots) {
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
        if (AnalysisConfig.cfgFillQAPlots) {
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
      auto posDauTrack = k0sCand.template posTrack_as<TrackCandidates>();
      auto negDauTrack = k0sCand.template negTrack_as<TrackCandidates>();

      /// Daughters
      // Positve pion
      auto trkppt = posDauTrack.pt();
      auto trkpy = posDauTrack.y();
      auto istrkphasTOF = posDauTrack.hasTOF();
      auto trkpNSigmaPiTPC = posDauTrack.tpcNSigmaPi();
      auto trkpNSigmaPiTOF = (istrkphasTOF) ? posDauTrack.tofNSigmaPi() : -999.;
      // Negative pion
      auto trknpt = negDauTrack.pt();
      auto trkny = negDauTrack.y();
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

      lResoSecondary = LorentzVectorSetXYZM(k0sCand.px(), k0sCand.py(), k0sCand.pz(), trkkMass);
      auto lPhiMinusPsiK0s = RecoDecay::constrainAngle(lResoSecondary.Phi() - lEPDet, 0.0, 2); // constrain angle to range 0, Pi
      //      auto v2K0s = std::cos(static_cast<float>(nmode) * lPhiMinusPsiK0s);

      float cosNPhiK0s = std::cos(static_cast<float>(nmode) * lResoSecondary.Phi());
      float sinNPhiK0s = std::sin(static_cast<float>(nmode) * lResoSecondary.Phi());

      auto v2K0s = cosNPhiK0s * collision.qvecRe()[lQvecDetInd] + sinNPhiK0s * collision.qvecIm()[lQvecDetInd];
      if constexpr (!IsMix) {
        if (AnalysisConfig.cfgFillQAPlots) {
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
          histos.fill(HIST("QA/before/hy_Secondary"), trkky);
          histos.fill(HIST("QA/before/hDCAtoPVSecondary"), trkkDCAtoPV);
          histos.fill(HIST("QA/before/hCPASecondary"), trkkCPA);
          histos.fill(HIST("QA/before/hPropTauSecondary"), trkkPropTau);
          histos.fill(HIST("QA/before/hInvmassSecondary"), trkkMass);

          histos.fill(HIST("QA/before/k0sv2vsinvmass"), lResoSecondary.M(), v2K0s);
        }
      }

      if (!SecondaryCuts.cfgByPassDauPIDSelection && !selectionPIDPion(posDauTrack))
        continue;
      if (!SecondaryCuts.cfgByPassDauPIDSelection && !selectionPIDPion(negDauTrack))
        continue;
      if (!SecondaryCuts.cfgByPassDauRapiditySelection && std::fabs(trkpy) > SecondaryCuts.cfgSecondaryDauRapidityMax)
        continue;
      if (!SecondaryCuts.cfgByPassDauRapiditySelection && std::fabs(trkny) > SecondaryCuts.cfgSecondaryDauRapidityMax)
        continue;
      if (!selectionK0s(collision, k0sCand))
        continue;

      if constexpr (!IsMix) {
        if (AnalysisConfig.cfgFillQAPlots) {
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
          histos.fill(HIST("QA/after/hy_Secondary"), trkky);
          histos.fill(HIST("QA/after/hDCAtoPVSecondary"), trkkDCAtoPV);
          histos.fill(HIST("QA/after/hCPASecondary"), trkkCPA);
          histos.fill(HIST("QA/after/hPropTauSecondary"), trkkPropTau);
          histos.fill(HIST("QA/after/hInvmassSecondary"), trkkMass);

          histos.fill(HIST("QA/after/k0sv2vsinvmass"), lResoSecondary.M(), v2K0s);
          if (AnalysisConfig.cfgFillAdditionalAxis) {
            if (AnalysisConfig.cfgFillOccupancy) {
              histos.fill(HIST("hInvmass_K0s"), lCentrality, lResoSecondary.Pt(), lResoSecondary.M(), v2K0s, static_cast<float>(nmode) * lPhiMinusPsiK0s, collision.trackOccupancyInTimeRange());
            } else {
              histos.fill(HIST("hInvmass_K0s"), lCentrality, lResoSecondary.Pt(), lResoSecondary.M(), v2K0s, static_cast<float>(nmode) * lPhiMinusPsiK0s);
            }
          } else {
            if (AnalysisConfig.cfgFillOccupancy) {
              histos.fill(HIST("hInvmass_K0s"), lCentrality, lResoSecondary.Pt(), lResoSecondary.M(), v2K0s, collision.trackOccupancyInTimeRange());
            } else {
              histos.fill(HIST("hInvmass_K0s"), lCentrality, lResoSecondary.Pt(), lResoSecondary.M(), v2K0s);
            }
          }
        }
        k0sIndicies.push_back(k0sCand.index());
      }
    }

    for (const auto& trackIndex : trackIndicies) {
      for (const auto& k0sIndex : k0sIndicies) {
        auto bTrack = dTracks1.rawIteratorAt(trackIndex);
        auto k0sCand = dTracks2.rawIteratorAt(k0sIndex);
        auto trkkMass = k0sCand.mK0Short();

        lDecayDaughter_bach = LorentzVectorSetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), MassPionCharged);
        lResoSecondary = LorentzVectorSetXYZM(k0sCand.px(), k0sCand.py(), k0sCand.pz(), trkkMass);
        lResoKstar = lResoSecondary + lDecayDaughter_bach;
        auto resoPhi = lResoKstar.Phi();
        // EP method
        auto lPhiMinusPsiKstar = RecoDecay::constrainAngle(resoPhi - lEPDet, 0.0, 2); // constrain angle to range 0, Pi
        auto resoFlowValue = std::cos(static_cast<float>(nmode) * lPhiMinusPsiKstar);
        // Scalar product method
        if (AnalysisConfig.cfgUseScalProduct) {
          float cosNPhi = std::cos(static_cast<float>(nmode) * resoPhi);
          float sinNPhi = std::sin(static_cast<float>(nmode) * resoPhi);
          resoFlowValue = cosNPhi * collision.qvecRe()[lQvecDetInd] + sinNPhi * collision.qvecIm()[lQvecDetInd];
        }

        // QA plots
        if constexpr (!IsMix) {
          if (AnalysisConfig.cfgFillQAPlots) {
            histos.fill(HIST("QA/before/KstarRapidity"), lResoKstar.Rapidity());
            histos.fill(HIST("QA/before/kstarinvmass"), lResoKstar.M());
            histos.fill(HIST("QA/before/kstarv2vsinvmass"), lResoKstar.M(), resoFlowValue);
          }
        }

        if (lResoKstar.Rapidity() > KstarCuts.cfgKstarMaxRap || lResoKstar.Rapidity() < KstarCuts.cfgKstarMinRap)
          continue;

        if constexpr (!IsMix) {
          unsigned int typeKstar = bTrack.sign() > 0 ? BinType::kKstarP : BinType::kKstarN;

          if (AnalysisConfig.cfgFillQAPlots) {
            histos.fill(HIST("QA/after/KstarRapidity"), lResoKstar.Rapidity());
            histos.fill(HIST("QA/after/kstarinvmass"), lResoKstar.M());
            histos.fill(HIST("QA/after/kstarv2vsinvmass"), lResoKstar.M(), resoFlowValue);
          }
          if (AnalysisConfig.cfgFillAdditionalAxis) {
            if (AnalysisConfig.cfgFillOccupancy) {
              histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResoKstar.Pt(), lResoKstar.M(), resoFlowValue, static_cast<float>(nmode) * lPhiMinusPsiKstar, collision.trackOccupancyInTimeRange());
            } else {
              histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResoKstar.Pt(), lResoKstar.M(), resoFlowValue, static_cast<float>(nmode) * lPhiMinusPsiKstar);
            }
          } else {
            if (AnalysisConfig.cfgFillOccupancy) {
              histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResoKstar.Pt(), lResoKstar.M(), resoFlowValue, collision.trackOccupancyInTimeRange());
            } else {
              histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResoKstar.Pt(), lResoKstar.M(), resoFlowValue);
            }
          }

          if (BkgEstimationConfig.cfgFillRotBkg) {
            for (int i = 0; i < BkgEstimationConfig.cfgNrotBkg; i++) {
              auto lRotAngle = BkgEstimationConfig.cfgMinRot + i * ((BkgEstimationConfig.cfgMaxRot - BkgEstimationConfig.cfgMinRot) / (BkgEstimationConfig.cfgNrotBkg - 1));
              if (AnalysisConfig.cfgFillQAPlots) {
                histos.fill(HIST("QA/RotBkg/hRotBkg"), lRotAngle);
              }
              if (BkgEstimationConfig.cfgRotPion) {
                lDaughterRot = lDecayDaughter_bach;
                ROOT::Math::RotationZ rot(lRotAngle);
                auto p3 = rot * lDaughterRot.Vect();
                lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                lResonanceRot = lDaughterRot + lResoSecondary;
              } else {
                lDaughterRot = lResoSecondary;
                ROOT::Math::RotationZ rot(lRotAngle);
                auto p3 = rot * lDaughterRot.Vect();
                lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                lResonanceRot = lDecayDaughter_bach + lDaughterRot;
              }
              resoPhi = lResonanceRot.Phi();
              auto lPhiMinusPsiKstar = RecoDecay::constrainAngle(resoPhi - lEPDet, 0.0, 2); // constrain angle to range 0, Pi
              auto resoFlowValue = std::cos(static_cast<float>(nmode) * lPhiMinusPsiKstar);
              if (AnalysisConfig.cfgUseScalProduct) {
                float cosNPhi = std::cos(static_cast<float>(nmode) * resoPhi);
                float sinNPhi = std::sin(static_cast<float>(nmode) * resoPhi);
                resoFlowValue = cosNPhi * collision.qvecRe()[lQvecDetInd] + sinNPhi * collision.qvecIm()[lQvecDetInd];
              }
              typeKstar = bTrack.sign() > 0 ? BinType::kKstarP_Rot : BinType::kKstarN_Rot;
              if (AnalysisConfig.cfgFillAdditionalAxis) {
                if (AnalysisConfig.cfgFillOccupancy) {
                  histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResonanceRot.Pt(), lResonanceRot.M(), resoFlowValue, static_cast<float>(nmode) * lPhiMinusPsiKstar, collision.trackOccupancyInTimeRange());
                } else {
                  histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResonanceRot.Pt(), lResonanceRot.M(), resoFlowValue, static_cast<float>(nmode) * lPhiMinusPsiKstar);
                }
              } else {
                if (AnalysisConfig.cfgFillOccupancy) {
                  histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResonanceRot.Pt(), lResonanceRot.M(), resoFlowValue, collision.trackOccupancyInTimeRange());
                } else {
                  histos.fill(HIST("hInvmass_Kstar"), typeKstar, lCentrality, lResonanceRot.Pt(), lResonanceRot.M(), resoFlowValue);
                }
              }
            }
          }
        } // IsMix

      } // k0sCand
    } // bTrack

    count++;

  } // fillHistograms

  // process dummy
  void processDummy(aod::Collisions const&)
  {
  }
  PROCESS_SWITCH(Chk892Flow, processDummy, "process Dummy", false);

  // process data
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
    float qAmpThr = 1e-4f;
    if (AnalysisConfig.cfgQvecSel && (collision.qvecAmp()[lDetId] < qAmpThr || collision.qvecAmp()[lRefAId] < qAmpThr || collision.qvecAmp()[lRefBId] < qAmpThr))
      return; // If we don't have a Q-vector
    lCentrality = getCentrality(collision);
    // lCentrality = collision.centFT0C();
    if (lCentrality < EventCuts.cfgEventCentralityMin || lCentrality > EventCuts.cfgEventCentralityMax)
      return;
    colCuts.fillQA(collision);

    fillHistograms<false, false, false>(collision, tracks, v0s, EventPlaneConfig.cfgnMods); // second order
  }
  PROCESS_SWITCH(Chk892Flow, processData, "Process Event for data without Partitioning", true);

  void processMC(aod::McCollisions const&)
  {
  }
  PROCESS_SWITCH(Chk892Flow, processMC, "Process MC for efficiency correction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Chk892Flow>(cfgc)};
}
