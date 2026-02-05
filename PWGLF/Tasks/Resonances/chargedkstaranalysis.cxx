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

/// \file chargedkstaranalysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
///
/// \author Navneet Kumar <navneet.kumar@cern.ch>, Protay

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/RCTSelectionFlags.h"
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
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector2.h"
// #include <TDatabasePDG.h> // FIXME
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h> // FIXME

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::aod::rctsel;

struct Chargedkstaranalysis {

  enum CentralityEstimator : int {
    FT0C = 1,
    FT0M = 2
  };
  struct : ConfigurableGroup {
    ConfigurableAxis cfgvtxbins{"cfgvtxbins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgmultbins{"cfgmultbins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};
    ConfigurableAxis cfgBinsPt{"cfgBinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsPtQA{"cfgBinsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsCent{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
    ConfigurableAxis cfgBinsVtxZ{"cfgBinsVtxZ", {VARIABLE_WIDTH, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "Binning of the z-vertex axis"};

    ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};
    ConfigurableAxis configThnAxisPhi{"configThnAxisPhi", {70, 0.0f, 7.0f}, "Phi axis"}; // 0 to 2pi

  } axisCfgs;

  struct : ConfigurableGroup {

    Configurable<bool> qAOptimisation{"qAOptimisation", false, "QA for optimisation with multiple THnSparse Axes"};
    // output THnSparses
    Configurable<bool> activateHelicityFrame{"activateHelicityFrame", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
    Configurable<bool> activateCollinsSoperFrame{"activateCollinsSoperFrame", false, "Activate the THnSparse with cosThStar w.r.t. Collins soper axis"};
    Configurable<bool> activateProductionFrame{"activateProductionFrame", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
    Configurable<bool> activateBeamAxisFrame{"activateBeamAxisFrame", false, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
    Configurable<bool> activateRandomFrame{"activateRandomFrame", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
    Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};
    Configurable<bool> cBoostKShot{"cBoostKShot", true, "Boost the Kshot in Charged Kstar frame of reference"};
    // Other cuts on Ks
    Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};

    // fixed variables
    float rapidityMotherData = 0.5;
    float beamEnergy = 13600.0;
    double beamMomentum = std::sqrt(beamEnergy * beamEnergy / 4 - o2::constants::physics::MassProton * o2::constants::physics::MassProton); // GeV
    int noOfDaughters = 2;
  } helicityCfgs;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  //  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::Mults>;
  //  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCPi, aod::pidTP     CKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  //    using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  using V0Candidates = aod::V0Datas;

  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Join<TrackCandidates, aod::McTrackLabels>;
  using MCV0Candidates = soa::Join<V0Candidates, aod::McV0Labels>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hChaKstar{"hChaKstar", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  // Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};

  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQ     ualityTracks, 5 -> kInAcceptanceTracks"};
  //
  //  Filter trackFilter = (trackSelection.node() == 0) || // from tpcSkimsTableCreator
  //                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
  //                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
  //                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
  //                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
  //                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));
  // Filter trackEtaFilter = nabs(aod::track::eta) < cfgCutEta; // Eta cut
  //
  // Configurables
  Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};

  struct RCTCut : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};

    RCTFlagsChecker rctChecker;
  };

  RCTCut rctCut;

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  struct : ConfigurableGroup {
    Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<int> confEvtOccupancyInTimeRangeMax{"confEvtOccupancyInTimeRangeMax", -1, "Evt sel: maximum track occupancy"};
    Configurable<int> confEvtOccupancyInTimeRangeMin{"confEvtOccupancyInTimeRangeMin", -1, "Evt sel: minimum track occupancy"};
    Configurable<bool> confEvtTriggerCheck{"confEvtTriggerCheck", false, "Evt sel: check for trigger"};
    Configurable<bool> confEvtOfflineCheck{"confEvtOfflineCheck", true, "Evt sel: check for offline selection"};
    Configurable<bool> confEvtTriggerTVXSel{"confEvtTriggerTVXSel", false, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> confEvtTFBorderCut{"confEvtTFBorderCut", false, "Evt sel: apply TF border cut"};
    Configurable<bool> confEvtUseITSTPCvertex{"confEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
    Configurable<bool> confEvtZvertexTimedifference{"confEvtZvertexTimedifference", true, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> confEvtPileupRejection{"confEvtPileupRejection", true, "Evt sel: apply pileup rejection"};
    Configurable<bool> confEvtNoITSROBorderCut{"confEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> confincludeCentralityMC{"confincludeCentralityMC", false, "Include centrality in MC"};
    Configurable<bool> confEvtCollInTimeRangeStandard{"confEvtCollInTimeRangeStandard", true, "Evt sel: apply NoCollInTimeRangeStandard"};
    Configurable<bool> confEvtCollInTimeRangeNarrow{"confEvtCollInTimeRangeNarrow", false, "Evt sel: apply NoCollInTimeRangeNarrow"};
    Configurable<bool> confEvtCollInTimeRangeStrict{"confEvtCollInTimeRangeStrict", false, "Evt sel: apply NoCollInTimeRangeStrict"};
    Configurable<bool> confEvtNoCollInRofStandard{"confEvtNoCollInRofStandard", false, "Evt sel: apply NoCollInRofStandard"};
    Configurable<bool> confEvtNoCollInRofStrict{"confEvtNoCollInRofStrict", false, "Evt sel: apply NoCollInRofStrict"};
    Configurable<bool> confEvtNoHighMultCollInPrevRof{"confEvtNoHighMultCollInPrevRof", false, "Evt sel: apply NoHighMultCollInPrevRof"};
    Configurable<bool> confEvtGoodITSLayersAll{"confEvtGoodITSLayersAll", false, "Evt sel: require all good ITS layers"};
    Configurable<bool> confEvtGoodITSLayer3{"confEvtGoodITSLayer3", false, "Evt sel: require good ITS layer 3"};
    Configurable<bool> confEvtGoodITSLayer0123{"confEvtGoodITSLayer0123", false, "Evt sel: require good ITS layers 0–3"};
    Configurable<bool> confEvtVertexTRDmatched{"confEvtVertexTRDmatched", false, "Evt sel: require TRD-matched vertex"};
    Configurable<bool> confEvtBBT0A{"confEvtBBT0A", false, "Evt sel: apply BBT0A cut"};
    Configurable<bool> confEvtBBT0C{"confEvtBBT0C", false, "Evt sel: apply BBT0C cut"};

    Configurable<float> cfgEventCentralityMin{"cfgEventCentralityMin", 0.0f, "Event sel: minimum centrality"};
    Configurable<float> cfgEventCentralityMax{"cfgEventCentralityMax", 100.0f, "Event sel: maximum centrality"};
    Configurable<int> cfgCentEst{"cfgCentEst", static_cast<int>(CentralityEstimator::FT0C), "Centrality estimator: 1=FT0C, 2=FT0M"};
  } eventCutCfgs;

  /// Track selections
  //
  struct : ConfigurableGroup {
    Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
    Configurable<float> cMaxEtacut{"cMaxEtacut", 0.8, "Track maximum eta cut"};
    // DCAr to PV
    Configurable<float> cMaxbDCArToPVcut{"cMaxbDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
    // DCAz to PV
    Configurable<float> cMaxbDCAzToPVcut{"cMaxbDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};

    // Track selections
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor

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
    Configurable<bool> cfgpTdepDCAxyCut{"cfgpTdepDCAxyCut", true, "pT-dependent DCAxy cut"};
    Configurable<bool> cfgpTdepDCAzCut{"cfgpTdepDCAzCut", false, "pT-dependent DCAz cut"};
    Configurable<float> cfgMaxbDCArToPVcut{"cfgMaxbDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
    Configurable<float> cfgMaxbDCAzToPVcut{"cfgMaxbDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
  } trackCutCfgs;

  struct : ConfigurableGroup {

    /// PID Selections, pion
    Configurable<bool> cfgTPConly{"cfgTPConly", false, "Use only TPC for PID"};                                     // bool
    Configurable<float> cfgMaxTPCnSigmaPion{"cfgMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};                 // TPC
    Configurable<float> cfgMaxTOFnSigmaPion{"cfgMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                 // TOF
    Configurable<float> cfgNsigmaCutCombinedPion{"cfgNsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"}; // Combined
    Configurable<bool> cfgTOFVeto{"cfgTOFVeto", false, "TOF Veto, if false, TOF is nessessary for PID selection"};  // TOF Veto
    Configurable<float> cfgTOFMinPt{"cfgTOFMinPt", 0.6, "Minimum TOF pT cut for Pion"};                             // TOF pT cut
  } pidCutCfgs;

  struct : ConfigurableGroup {

    // Secondary Selection
    Configurable<bool> cfgByPassDauPIDSelection{"cfgByPassDauPIDSelection", true, "Bypass Daughters PID selection"};
    Configurable<float> cSecondaryDauDCAMax{"cSecondaryDauDCAMax", 1., "Maximum DCA Secondary daughters to PV"};
    Configurable<float> cSecondaryDauPosDCAtoPVMin{"cSecondaryDauPosDCAtoPVMin", 0.0, "Minimum DCA Secondary positive daughters to PV"};
    Configurable<float> cSecondaryDauNegDCAtoPVMin{"cSecondaryDauNegDCAtoPVMin", 0.0, "Minimum DCA Secondary negative daughters to PV"};

    Configurable<float> cSecondaryPtMin{"cSecondaryPtMin", 0.f, "Minimum transverse momentum of Secondary"};
    Configurable<float> cSecondaryRapidityMax{"cSecondaryRapidityMax", 0.5, "Maximum rapidity of Secondary"};
    Configurable<float> cSecondaryRadiusMin{"cSecondaryRadiusMin", 0, "Minimum transverse radius of Secondary"};
    Configurable<float> cSecondaryRadiusMax{"cSecondaryRadiusMax", 999.9, "Maximum transverse radius of Secondary"};
    Configurable<float> cSecondaryCosPAMin{"cSecondaryCosPAMin", 0.995, "Mininum cosine pointing angle of Secondary"};
    Configurable<float> cSecondaryDCAtoPVMax{"cSecondaryDCAtoPVMax", 0.3, "Maximum DCA Secondary to PV"};
    Configurable<float> cSecondaryProperLifetimeMax{"cSecondaryProperLifetimeMax", 20, "Maximum Secondary Lifetime"};
    Configurable<float> cSecondaryMassWindow{"cSecondaryMassWindow", 0.075, "Secondary inv mass selciton window"};
    Configurable<float> cfgSecondaryparamArmenterosCut{"cfgSecondaryparamArmenterosCut", 0.2, "parameter for Armenteros Cut"};
    Configurable<bool> cfgSecondaryCrossMassHypothesisCut{"cfgSecondaryCrossMassHypothesisCut", false, "Apply cut based on the lambda mass hypothesis"};
    Configurable<float> cfgSecondaryCrossMassCutWindow{"cfgSecondaryCrossMassCutWindow", 0.05, "Secondary inv mass selection window with (anti)lambda hypothesis"};

  } secondaryCutsCfgs;

  // K* selection
  struct : ConfigurableGroup {
    Configurable<float> cKstarMaxRap{"cKstarMaxRap", 0.5, "Kstar maximum rapidity"};
    Configurable<float> cKstarMinRap{"cKstarMinRap", -0.5, "Kstar minimum rapidity"};
  } kstarCutCfgs;

  struct : ConfigurableGroup {
    // For rotational background
    Configurable<bool> cFillRotBkg{"cFillRotBkg", true, "Fill rotated background"};
    Configurable<float> confMinRot{"confMinRot", 5.0 * o2::constants::math::PI / 6.0, "Minimum of rotation"};
    Configurable<float> confMaxRot{"confMaxRot", 7.0 * o2::constants::math::PI / 6.0, "Maximum of rotation"};
    Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  } rotBkgEstCfgs;

  float centrality;

  // PDG code
  int kPDGK0s = kK0Short;
  int kKstarPlus = static_cast<int>(o2::constants::physics::Pdg::kKPlusStar892);
  // int kPiPlus = kPiPlus;
  int kPDGK0 = kK0;
  // Variable declaration
  ROOT::Math::PxPyPzEVector beam1{0., 0., -helicityCfgs.beamMomentum, 13600. / 2.};
  ROOT::Math::PxPyPzEVector beam2{0., 0., helicityCfgs.beamMomentum, 13600. / 2.};

  void init(o2::framework::InitContext&)
  {
    centrality = -999;

    rctCut.rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    colCuts.setCuts(eventCutCfgs.confEvtZvtx, eventCutCfgs.confEvtTriggerCheck, eventCutCfgs.confEvtOfflineCheck, /*checkRun3*/ true, /*triggerTVXsel*/ false, eventCutCfgs.confEvtOccupancyInTimeRangeMax, eventCutCfgs.confEvtOccupancyInTimeRangeMin);
    colCuts.init(&histos);
    colCuts.setTriggerTVX(eventCutCfgs.confEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(eventCutCfgs.confEvtTFBorderCut);
    colCuts.setApplyNoITSROBorderCut(eventCutCfgs.confEvtNoITSROBorderCut);
    colCuts.setApplyITSTPCvertex(eventCutCfgs.confEvtUseITSTPCvertex);
    colCuts.setApplyZvertexTimedifference(eventCutCfgs.confEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(eventCutCfgs.confEvtPileupRejection);
    colCuts.setApplyCollInTimeRangeStandard(eventCutCfgs.confEvtCollInTimeRangeStandard);
    colCuts.setApplyCollInTimeRangeNarrow(eventCutCfgs.confEvtCollInTimeRangeNarrow);
    colCuts.setApplyCollInTimeRangeStrict(eventCutCfgs.confEvtCollInTimeRangeStrict);
    colCuts.setApplyCollInRofStandard(eventCutCfgs.confEvtNoCollInRofStandard);
    colCuts.setApplyCollInRofStrict(eventCutCfgs.confEvtNoCollInRofStrict);
    colCuts.setApplyHighMultCollInPrevRof(eventCutCfgs.confEvtNoHighMultCollInPrevRof);
    colCuts.setApplyGoodITSLayersAll(eventCutCfgs.confEvtGoodITSLayersAll);
    colCuts.setApplyGoodITSLayer3(eventCutCfgs.confEvtGoodITSLayer3);
    colCuts.setApplyGoodITSLayer0123(eventCutCfgs.confEvtGoodITSLayer0123);
    colCuts.setApplyVertexTRDmatched(eventCutCfgs.confEvtVertexTRDmatched);
    colCuts.setApplyBBT0A(eventCutCfgs.confEvtBBT0A);
    colCuts.setApplyBBT0C(eventCutCfgs.confEvtBBT0C);

    AxisSpec centAxis = {axisCfgs.cfgBinsCent, "T0M (%)"};
    AxisSpec vtxzAxis = {axisCfgs.cfgBinsVtxZ, "Z Vertex (cm)"};
    AxisSpec ptAxis = {axisCfgs.cfgBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {axisCfgs.cfgBinsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec radiusAxis = {50, 0, 5, "Radius (cm)"};
    AxisSpec cpaAxis = {50, 0.95, 1.0, "CPA"};
    AxisSpec tauAxis = {250, 0, 25, "Lifetime (cm)"};
    AxisSpec dcaAxis = {200, 0, 2, "DCA (cm)"};
    AxisSpec dcaxyAxis = {200, 0, 2, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {200, 0, 2, "DCA_{#it{z}} (cm)"};
    AxisSpec yAxis = {100, -1, 1, "Rapidity"};
    AxisSpec invMassAxisK0s = {400 / cNbinsDiv, 0.3, 0.7, "Invariant Mass (GeV/#it{c}^2)"};    // K0s ~497.611
    AxisSpec invMassAxisReso = {900 / cNbinsDiv, 0.5f, 1.4f, "Invariant Mass (GeV/#it{c}^2)"}; // chK(892) ~892
    AxisSpec invMassAxisScan = {150, 0, 1.5, "Invariant Mass (GeV/#it{c}^2)"};                 // For selection
    AxisSpec pidQAAxis = {130, -6.5, 6.5};
    AxisSpec dataTypeAxis = {9, 0, 9, "Histogram types"};
    AxisSpec mcTypeAxis = {4, 0, 4, "Histogram types"};
    AxisSpec thnAxisPOL{axisCfgs.configThnAxisPOL, "Configurabel theta axis"};
    AxisSpec thnAxisPhi = {axisCfgs.configThnAxisPhi, "Configurabel phi axis"}; // 0 to 2pi
    // THnSparse
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};

    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{5, 0, 5.0}});
    std::shared_ptr<TH1> hCutFlow = histos.get<TH1>(HIST("hEvtSelInfo"));
    hCutFlow->GetXaxis()->SetBinLabel(1, "All Events");
    hCutFlow->GetXaxis()->SetBinLabel(2, "coll cuts");
    hCutFlow->GetXaxis()->SetBinLabel(3, "rctChecker");
    hCutFlow->GetXaxis()->SetBinLabel(4, "Multiplicity");
    hCutFlow->GetXaxis()->SetBinLabel(5, "IsINELgt0");

    constexpr int kNTrackCuts = 22;

    histos.add("QA/hTrackCutFlow", "Track cut flow", kTH1I, {{kNTrackCuts, 0.5, kNTrackCuts + 0.5}});

    auto hTrackCutFlow = histos.get<TH1>(HIST("QA/hTrackCutFlow"));

    int bin = 1;
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "All tracks");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "pT min");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "|eta| max");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "ITS clusters");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "TPC clusters");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "TPC crossed rows ratio");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "ITS chi2/Ncl");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "TPC chi2/Ncl");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "Has ITS");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "Has TPC");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "Has TOF");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "ITS refit");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "TPC refit");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "PV contributor");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "Global w/o DCA");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "Global track");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "Primary track");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "DCAxy max");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "DCAz max");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "pT-dep DCAxy");
    hTrackCutFlow->GetXaxis()->SetBinLabel(bin++, "pT-dep DCAz");

    constexpr int kNK0sCuts = 14;
    int iK0sbin = 1;
    histos.add("QA/K0sCutCheck", "K0s cut flow", kTH1I, {{kNK0sCuts, 0.5, kNK0sCuts + 0.5}});
    auto hK0sCut = histos.get<TH1>(HIST("QA/K0sCutCheck"));
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "All PASS");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "DauDCA>max");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "PosDCAtoPV<min");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "NegDCAtoPV<min");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "pT<min");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "|y|>max");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "R<min or R>max");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "DCAtoPV>max");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "cosPA<min");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "ctau>max");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "qtarm<a|alpha|");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "|M(K0s)-m0|>win");
    hK0sCut->GetXaxis()->SetBinLabel(iK0sbin++, "cross-mass veto");

    histos.add("QA/before/CentDist", "Centrality distribution", {HistType::kTH1D, {centAxis}});
    histos.add("QA/before/CentDist1", "Centrality distribution", o2::framework::kTH1F, {{110, 0, 110}});
    histos.add("QA/before/VtxZ", "Centrality distribution", {HistType::kTH1D, {vtxzAxis}});
    histos.add("QA/before/hEvent", "Number of Events", HistType::kTH1F, {{1, 0.5, 1.5}});

    histos.add("QA/trkbpionTPCPIDME", "TPC PID of bachelor pion candidates", HistType::kTH2D, {ptAxisQA, pidQAAxis});

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
    histos.add("hInvmass_Kstar", "Invariant mass of unlike-sign chK(892)", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_KstarME", "Invariant mass of unlike-sign chK(892)ME", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_KstarRotated", "Invariant mass of unlike-sign chK(892)Rota", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});

    // Mass QA (quick check)
    histos.add("QA/before/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
    histos.add("QA/before/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

    histos.add("QA/after/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
    histos.add("QA/after/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
    histos.add("QA/after/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

    if (!helicityCfgs.qAOptimisation) {
      hChaKstar.add("h3ChaKstarInvMassDS", "h3ChaKstarInvMassDS", kTHnSparseF, {centAxis, ptAxis, invMassAxisReso, thnAxisPOL, thnAxisPhi}, true);
      hChaKstar.add("h3ChaKstarInvMassME", "h3ChaKstarInvMassME", kTHnSparseF, {centAxis, ptAxis, invMassAxisReso, thnAxisPOL, thnAxisPhi}, true);
      hChaKstar.add("h3ChaKstarInvMassRot", "h3ChaKstarInvMassRot", kTHnSparseF, {centAxis, ptAxis, invMassAxisReso, thnAxisPOL, thnAxisPhi}, true);
    }

    if (rotBkgEstCfgs.cFillRotBkg) {
      histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});
    }
    // MC
    if (doprocessMC) {

      histos.add("QAMC/hEvent", "Number of Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      // Bachelor pion
      histos.add("QAMC/trkbpionDCAxy", "DCAxy distribution of bachelor pion candidates", HistType::kTH1D, {dcaxyAxis});
      histos.add("QAMC/trkbpionDCAz", "DCAz distribution of bachelor pion candidates", HistType::kTH1D, {dcazAxis});
      histos.add("QAMC/trkbpionpT", "pT distribution of bachelor pion candidates", HistType::kTH1D, {ptAxis});
      histos.add("QAMC/trkbpionTPCPID", "TPC PID of bachelor pion candidates", HistType::kTH2D, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkbpionTOFPID", "TOF PID of bachelor pion candidates", HistType::kTH2D, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkbpionTPCTOFPID", "TPC-TOF PID map of bachelor pion candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});

      // Secondary pion 1
      histos.add("QAMC/trkppionDCAxy", "DCAxy distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcaxyAxis});
      histos.add("QAMC/trkppionDCAz", "DCAz distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {dcazAxis});
      histos.add("QAMC/trkppionpT", "pT distribution of secondary pion 1 (positive) candidates", HistType::kTH1D, {ptAxis});
      histos.add("QAMC/trkppionTPCPID", "TPC PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTOFPID", "TOF PID of secondary pion 1 (positive) candidates", HistType::kTH2D, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTPCTOFPID", "TPC-TOF PID map of secondary pion 1 (positive) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});

      // Secondary pion 2
      histos.add("QAMC/trknpionTPCPID", "TPC PID of secondary pion 2 (negative) candidates", HistType::kTH2D, {ptAxis, pidQAAxis});
      histos.add("QAMC/trknpionTOFPID", "TOF PID of secondary pion 2 (negative) candidates", HistType::kTH2D, {ptAxis, pidQAAxis});
      histos.add("QAMC/trknpionTPCTOFPID", "TPC-TOF PID map of secondary pion 2 (negative) candidates", HistType::kTH2D, {pidQAAxis, pidQAAxis});
      histos.add("QAMC/trknpionpT", "pT distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {ptAxis});
      histos.add("QAMC/trknpionDCAxy", "DCAxy distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {dcaxyAxis});
      histos.add("QAMC/trknpionDCAz", "DCAz distribution of secondary pion 2 (negative) candidates", HistType::kTH1D, {dcazAxis});

      // Secondary Resonance (K0s cand)
      histos.add("QAMC/hDauDCASecondary", "DCA of daughters of secondary resonance", HistType::kTH1D, {dcaAxis});
      histos.add("QAMC/hDauPosDCAtoPVSecondary", "Pos DCA to PV of daughters secondary resonance", HistType::kTH1D, {dcaAxis});
      histos.add("QAMC/hDauNegDCAtoPVSecondary", "Neg DCA to PV of daughters secondary resonance", HistType::kTH1D, {dcaAxis});

      histos.add("QAMC/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1D, {ptAxis});
      histos.add("QAMC/hy_Secondary", "Rapidity distribution of secondary resonance", HistType::kTH1D, {yAxis});
      histos.add("QAMC/hRadiusSecondary", "Radius distribution of secondary resonance", HistType::kTH1D, {radiusAxis});
      histos.add("QAMC/hCPASecondary", "Cosine pointing angle distribution of secondary resonance", HistType::kTH1D, {cpaAxis});
      histos.add("QAMC/hDCAtoPVSecondary", "DCA to PV distribution of secondary resonance", HistType::kTH1D, {dcaAxis});
      histos.add("QAMC/hPropTauSecondary", "Proper Lifetime distribution of secondary resonance", HistType::kTH1D, {tauAxis});
      histos.add("QAMC/hPtAsymSecondary", "pT asymmetry distribution of secondary resonance", HistType::kTH1D, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
      histos.add("QAMC/hInvmassSecondary", "Invariant mass of unlike-sign secondary resonance", HistType::kTH1D, {invMassAxisK0s});

      // K892
      histos.add("QAMC/KstarOA", "Opening angle of chK(892)", HistType::kTH1D, {AxisSpec{100, 0, 3.14, "Opening angle"}});
      histos.add("QAMC/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});

      histos.add("QAMC/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
      histos.add("QAMC/kstarinvmass_noKstar", "Invariant mass of unlike-sign no chK(892)", HistType::kTH1D, {invMassAxisReso});

      histos.add("hInvmass_Kstar_MC", "Invariant mass of unlike chK(892)", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});

      ccdb->setURL(cfgURL);
      ccdbApi.init("http://alice-ccdb.cern.ch");
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    }

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in chK(892) Analysis Task";
    histos.print();
  }
  float lMultiplicity;
  template <typename CollisionType>
  float getCentrality(CollisionType const& collision)
  {
    auto centEst = static_cast<CentralityEstimator>(eventCutCfgs.cfgCentEst.value);

    switch (centEst) {
      case CentralityEstimator::FT0C:
        return collision.centFT0C();
      case CentralityEstimator::FT0M:
        return collision.centFT0M();
      default:
        return -999.f;
    }
  }

  // Track selection
  template <typename TrackType>
  bool trackCut(TrackType const& track)
  {
    int ibin = 1;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    // basic track cuts
    if (std::abs(track.pt()) < trackCutCfgs.cMinPtcut)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (std::abs(track.eta()) > trackCutCfgs.cMaxEtacut)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (track.itsNCls() < trackCutCfgs.cfgITScluster)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (track.tpcNClsFound() < trackCutCfgs.cfgTPCcluster)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (track.tpcCrossedRowsOverFindableCls() < trackCutCfgs.cfgRatioTPCRowsOverFindableCls)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (track.itsChi2NCl() >= trackCutCfgs.cfgITSChi2NCl)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (track.tpcChi2NCl() >= trackCutCfgs.cfgTPCChi2NCl)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgHasITS && !track.hasITS())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgHasTPC && !track.hasTPC())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgHasTOF && !track.hasTOF())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgUseITSRefit && !track.passedITSRefit())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgPVContributor && !track.isPVContributor())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (std::abs(track.dcaXY()) > trackCutCfgs.cMaxbDCArToPVcut)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (std::abs(track.dcaZ()) > trackCutCfgs.cMaxbDCAzToPVcut)
      return false;
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgpTdepDCAxyCut) {
      if (std::abs(track.dcaXY()) > (0.004 + (0.013 / track.pt())))
        return false;
    } else {
      if (std::abs(track.dcaXY()) > trackCutCfgs.cfgMaxbDCArToPVcut)
        return false;
    }
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);

    if (trackCutCfgs.cfgpTdepDCAzCut) {
      // Tuned on the LHC22f anchored MC LHC23d1d on primary pions. 7 Sigmas of the resolution
      if (std::abs(track.dcaZ()) > (0.004 + (0.013 / track.pt())))
        return false;
    } else {
      if (std::abs(track.dcaZ()) > trackCutCfgs.cfgMaxbDCAzToPVcut)
        return false;
    }
    histos.fill(HIST("QA/hTrackCutFlow"), ibin++);
    return true;
  }

  // PID selection tools
  template <typename TrackType>
  bool selectionPIDPion(TrackType const& candidate)
  {
    if (std::abs(candidate.tpcNSigmaPi()) >= pidCutCfgs.cfgMaxTPCnSigmaPion)
      return false;
    if (pidCutCfgs.cfgTPConly)
      return true;
    if (candidate.pt() <= pidCutCfgs.cfgTOFMinPt)
      return true;

    if (candidate.hasTOF()) {
      const bool tofPIDPassed = std::abs(candidate.tofNSigmaPi()) < pidCutCfgs.cfgMaxTOFnSigmaPion;
      const bool combo = (pidCutCfgs.cfgNsigmaCutCombinedPion > 0) &&
                         (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() +
                            candidate.tofNSigmaPi() * candidate.tofNSigmaPi() <
                          pidCutCfgs.cfgNsigmaCutCombinedPion * pidCutCfgs.cfgNsigmaCutCombinedPion);
      return tofPIDPassed || combo;
    } else {
      return pidCutCfgs.cfgTOFVeto;
    }
  }

  template <typename CollisionType, typename K0sType>
  bool selectionK0s(CollisionType const& collision, K0sType const& candidate)
  {
    int ibin = 1;
    auto dauDCA = candidate.dcaV0daughters();
    auto dauPosDCAtoPV = candidate.dcapostopv();
    auto dauNegDCAtoPV = candidate.dcanegtopv();
    auto pT = candidate.pt();
    auto rapidity = candidate.yK0Short();
    auto v0Radius = candidate.v0radius();
    auto dcaToPV = candidate.dcav0topv();
    auto cosPA = candidate.v0cosPA();
    auto propTauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0s;
    auto mK0s = candidate.mK0Short();
    auto mLambda = candidate.mLambda();
    auto mALambda = candidate.mAntiLambda();

    bool returnFlag = true;
    histos.fill(HIST("QA/K0sCutCheck"), ibin);
    if (dauDCA > secondaryCutsCfgs.cSecondaryDauDCAMax) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);
    if (dauPosDCAtoPV < secondaryCutsCfgs.cSecondaryDauPosDCAtoPVMin) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (dauNegDCAtoPV < secondaryCutsCfgs.cSecondaryDauNegDCAtoPVMin) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (pT < secondaryCutsCfgs.cSecondaryPtMin) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (std::fabs(rapidity) > secondaryCutsCfgs.cSecondaryRapidityMax) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (v0Radius < secondaryCutsCfgs.cSecondaryRadiusMin || v0Radius > secondaryCutsCfgs.cSecondaryRadiusMax) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (dcaToPV > secondaryCutsCfgs.cSecondaryDCAtoPVMax) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (cosPA < secondaryCutsCfgs.cSecondaryCosPAMin) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (propTauK0s > secondaryCutsCfgs.cSecondaryProperLifetimeMax) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (candidate.qtarm() < secondaryCutsCfgs.cfgSecondaryparamArmenterosCut * std::fabs(candidate.alpha())) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (std::fabs(mK0s - MassK0Short) > secondaryCutsCfgs.cSecondaryMassWindow) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);

    if (secondaryCutsCfgs.cfgSecondaryCrossMassHypothesisCut &&
        ((std::fabs(mLambda - MassLambda0) < secondaryCutsCfgs.cfgSecondaryCrossMassCutWindow) || (std::fabs(mALambda - MassLambda0Bar) < secondaryCutsCfgs.cfgSecondaryCrossMassCutWindow))) {
      returnFlag = false;
    }
    ibin++;
    if (returnFlag == true)
      histos.fill(HIST("QA/K0sCutCheck"), ibin);
    return returnFlag;

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

    if (std::abs(motherkV0.pdgCode()) != kPDGK0) // Is it K0s?
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
  double massPi = o2::constants::physics::MassPionCharged;
  double massK0s = o2::constants::physics::MassK0Short;

  template <typename T>
  void fillInvMass(const T& mother, float multiplicity, const T& daughter1, const T& daughter2, bool isMix)
  {
    TRandom* rn = new TRandom();
    rn->SetSeed(0);
    // Variable declarations
    ROOT::Math::PxPyPzMVector fourVecDauCM, daughterRot, motherRot, daughterRotCM;
    ROOT::Math::XYZVectorF beam1CM, beam2CM, zAxisCS, yAxisCS, xAxisCS;
    ROOT::Math::XYZVectorF v1CM, zaxisHE, yaxisHE, xaxisHE;
    ROOT::Math::XYZVector randomVec, beamVec, normalVec;
    float theta2;
    // //polarization calculations
    // zBeam = ROOT::Math::XYZVector(0.f, 0.f, 1.f); // ẑ: beam direction in lab frame

    ROOT::Math::Boost boost{mother.BoostToCM()}; // define the boost to the center of mass frame
    fourVecDauCM = boost(daughter1);             // boost the frame of daughter to the center of mass frame
    // threeVecDauCM = fourVecDauCM.Vect();         // get the 3 vector of daughter in the frame of mother

    beam1CM = ROOT::Math::XYZVectorF((boost(beam1).Vect()).Unit());
    beam2CM = ROOT::Math::XYZVectorF((boost(beam2).Vect()).Unit());

    v1CM = ROOT::Math::XYZVectorF(boost(daughter1).Vect()).Unit();
    // ROOT::Math::XYZVectorF v2_CM{(boost(daughter1).Vect()).Unit()};
    // using positive sign convention for the first track
    // ROOT::Math::XYZVectorF v_CM = (t1.sign() > 0 ? v1CM : v2_CM); // here selected decay daughter momentum is intested. here you can choose one decay daughter
    // Helicity Frame
    zaxisHE = ROOT::Math::XYZVectorF(mother.Vect()).Unit();
    yaxisHE = ROOT::Math::XYZVectorF(beam1CM.Cross(beam2CM)).Unit();
    xaxisHE = ROOT::Math::XYZVectorF(yaxisHE.Cross(zaxisHE)).Unit();

    // CosThetaHE = zaxisHE.Dot(v_CM);

    auto anglePhi = std::atan2(yaxisHE.Dot(v1CM), xaxisHE.Dot(v1CM));
    anglePhi = RecoDecay::constrainAngle(anglePhi, 0.0);
    // if (anglePhi < 0) {
    //   anglePhi += o2::constants::math::TwoPI; // ensure phi is in [0, 2pi]
    // }

    // CS Frame
    zAxisCS = ROOT::Math::XYZVectorF((beam1CM.Unit() - beam2CM.Unit())).Unit();
    yAxisCS = ROOT::Math::XYZVectorF(beam1CM.Cross(beam2CM)).Unit();
    xAxisCS = ROOT::Math::XYZVectorF(yAxisCS.Cross(zAxisCS)).Unit();
    double cosThetaStarCS = zAxisCS.Dot(v1CM);
    auto phiCS = std::atan2(yAxisCS.Dot(v1CM), xAxisCS.Dot(v1CM));
    phiCS = RecoDecay::constrainAngle(phiCS, 0.0);

    // if (std::abs(mother.Rapidity()) < config.rapidityMotherData) {
    if (helicityCfgs.activateHelicityFrame) {
      // helicityVec = mother.Vect(); // 3 vector of mother in COM frame
      // auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));
      auto cosThetaStarHelicity = mother.Vect().Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(mother.Vect().Mag2()));
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity, anglePhi);
        }

        for (int i = 0; i < helicityCfgs.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / helicityCfgs.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / helicityCfgs.rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

          motherRot = daughterRot + daughter2;

          ROOT::Math::Boost boost2{motherRot.BoostToCM()};
          daughterRotCM = boost2(daughterRot);

          auto cosThetaStarHelicityRot = motherRot.Vect().Dot(daughterRotCM.Vect()) / (std::sqrt(daughterRotCM.Vect().Mag2()) * std::sqrt(motherRot.Vect().Mag2()));
          auto phiHelicityRot = std::atan2(yaxisHE.Dot(daughterRotCM.Vect().Unit()), xaxisHE.Dot(daughterRotCM.Vect().Unit()));
          phiHelicityRot = RecoDecay::constrainAngle(phiHelicityRot, 0.0);
          if (motherRot.Rapidity() < helicityCfgs.rapidityMotherData)
            hChaKstar.fill(HIST("h3ChaKstarInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarHelicityRot, phiHelicityRot);
        }
      } else {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity, anglePhi);
        }
      }
    } else if (helicityCfgs.activateCollinsSoperFrame) {
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarCS, phiCS);
        }

        for (int i = 0; i < helicityCfgs.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / helicityCfgs.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / helicityCfgs.rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

          motherRot = daughterRot + daughter2;

          ROOT::Math::Boost boost2{motherRot.BoostToCM()};
          daughterRotCM = boost2(daughterRot);

          auto cosThetaStarCSrot = zAxisCS.Dot(daughterRotCM.Vect()) / std::sqrt(daughterRotCM.Vect().Mag2());
          auto phiCSrot = std::atan2(yAxisCS.Dot(daughterRotCM.Vect().Unit()), xAxisCS.Dot(daughterRotCM.Vect().Unit()));
          phiCSrot = RecoDecay::constrainAngle(phiCSrot, 0.0);

          if (motherRot.Rapidity() < helicityCfgs.rapidityMotherData)
            hChaKstar.fill(HIST("h3ChaKstarInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarCSrot, phiCSrot);
        }
      } else {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarCS, phiCS);
        }
      }
    } else if (helicityCfgs.activateProductionFrame) {
      normalVec = ROOT::Math::XYZVector(mother.Py(), -mother.Px(), 0.f);
      auto cosThetaProduction = normalVec.Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(normalVec.Mag2()));
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaProduction, anglePhi);
        }
        for (int i = 0; i < helicityCfgs.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / helicityCfgs.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / helicityCfgs.rotationalCut);
          motherRot = ROOT::Math::PxPyPzMVector(mother.Px() * std::cos(theta2) - mother.Py() * std::sin(theta2), mother.Px() * std::sin(theta2) + mother.Py() * std::cos(theta2), mother.Pz(), mother.M());
          if (std::abs(motherRot.Rapidity()) < helicityCfgs.rapidityMotherData) {
            hChaKstar.fill(HIST("h3ChaKstarInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaProduction, anglePhi);
          }
        }
      } else {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaProduction, anglePhi);
        }
      }
    } else if (helicityCfgs.activateBeamAxisFrame) {
      beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
      auto cosThetaStarBeam = beamVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam, anglePhi);
        }
        for (int i = 0; i < helicityCfgs.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / helicityCfgs.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / helicityCfgs.rotationalCut);
          motherRot = ROOT::Math::PxPyPzMVector(mother.Px() * std::cos(theta2) - mother.Py() * std::sin(theta2), mother.Px() * std::sin(theta2) + mother.Py() * std::cos(theta2), mother.Pz(), mother.M());
          if (std::abs(motherRot.Rapidity()) < helicityCfgs.rapidityMotherData) {
            hChaKstar.fill(HIST("h3ChaKstarInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarBeam, anglePhi);
          }
        }
      } else {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam, anglePhi);
        }
      }
    } else if (helicityCfgs.activateRandomFrame) {
      auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
      auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

      randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
      auto cosThetaStarRandom = randomVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());
      if (!isMix) {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassDS"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom, phiRandom);
        }
        for (int i = 0; i < helicityCfgs.cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / helicityCfgs.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / helicityCfgs.rotationalCut);
          motherRot = ROOT::Math::PxPyPzMVector(mother.Px() * std::cos(theta2) - mother.Py() * std::sin(theta2), mother.Px() * std::sin(theta2) + mother.Py() * std::cos(theta2), mother.Pz(), mother.M());
          if (std::abs(motherRot.Rapidity()) < helicityCfgs.rapidityMotherData) {
            hChaKstar.fill(HIST("h3ChaKstarInvMassRot"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarRandom, phiRandom);
          }
        }
      } else {
        if (std::abs(mother.Rapidity()) < helicityCfgs.rapidityMotherData) {
          hChaKstar.fill(HIST("h3ChaKstarInvMassME"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom, phiRandom);
        }
      }
    }
    // }
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType, typename TracksTypeK0s>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksTypeK0s& dTracks2)
  {
    histos.fill(HIST("QA/before/CentDist"), collision.centFT0M());
    histos.fill(HIST("QA/before/CentDist1"), collision.centFT0M());
    ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResoSecondary, lDecayDaughter_bach, lResoKstar, chargeKstarrot;
    std::vector<int> trackIndicies = {};
    std::vector<int> k0sIndicies = {};

    for (const auto& bTrack : dTracks1) {
      auto trkbpt = bTrack.pt();
      auto istrkbhasTOF = bTrack.hasTOF();
      auto trkbNSigmaPiTPC = bTrack.tpcNSigmaPi();
      auto trkbNSigmaPiTOF = (istrkbhasTOF) ? bTrack.tofNSigmaPi() : -999.;

      if constexpr (!IsMix) {
        // Bachelor pion QA plots
        histos.fill(HIST("QA/before/trkbpionTPCPID"), trkbpt, trkbNSigmaPiTPC);
        if (istrkbhasTOF) {
          histos.fill(HIST("QA/before/trkbpionTOFPID"), trkbpt, trkbNSigmaPiTOF);
          histos.fill(HIST("QA/before/trkbpionTPCTOFPID"), trkbNSigmaPiTPC, trkbNSigmaPiTOF);
        }
        histos.fill(HIST("QA/before/trkbpionpT"), trkbpt);
        histos.fill(HIST("QA/before/trkbpionDCAxy"), bTrack.dcaXY());
        histos.fill(HIST("QA/before/trkbpionDCAz"), bTrack.dcaZ());
      } else {

        histos.fill(HIST("QA/trkbpionTPCPIDME"), trkbpt, trkbNSigmaPiTPC);
      }

      if (!trackCut(bTrack))
        continue;
      if (!selectionPIDPion(bTrack))
        continue;

      if constexpr (!IsMix) {
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
      trackIndicies.push_back(bTrack.index());
    }

    for (const auto& K0scand : dTracks2) {
      auto posDauTrack = K0scand.template posTrack_as<TrackCandidates>();
      auto negDauTrack = K0scand.template negTrack_as<TrackCandidates>();

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
      auto trkkDauDCA = K0scand.dcaV0daughters();
      auto trkkDauDCAPostoPV = K0scand.dcapostopv();
      auto trkkDauDCANegtoPV = K0scand.dcanegtopv();
      auto trkkpt = K0scand.pt();
      auto trkky = K0scand.yK0Short();
      auto trkkRadius = K0scand.v0radius();
      auto trkkDCAtoPV = K0scand.dcav0topv();
      auto trkkCPA = K0scand.v0cosPA();
      auto trkkMass = K0scand.mK0Short();

      if constexpr (!IsMix) {
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
        histos.fill(HIST("QA/before/hInvmassSecondary"), trkkMass);
      }

      if (!secondaryCutsCfgs.cfgByPassDauPIDSelection && !selectionPIDPion(posDauTrack)) // Perhaps it's already applied in trackCut (need to check QA plots)
        continue;
      if (!secondaryCutsCfgs.cfgByPassDauPIDSelection && !selectionPIDPion(negDauTrack))
        continue;
      if (!selectionK0s(collision, K0scand))
        continue;

      if constexpr (!IsMix) {
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
        histos.fill(HIST("QA/after/hInvmassSecondary"), trkkMass);
      }
      k0sIndicies.push_back(K0scand.index());
    }

    for (const auto& trackIndex : trackIndicies) {
      for (const auto& k0sIndex : k0sIndicies) {
        auto bTrack = dTracks1.rawIteratorAt(trackIndex);
        auto k0Scand = dTracks2.rawIteratorAt(k0sIndex);
        lDecayDaughter_bach = ROOT::Math::PxPyPzMVector(bTrack.px(), bTrack.py(), bTrack.pz(), massPi);
        lResoSecondary = ROOT::Math::PxPyPzMVector(k0Scand.px(), k0Scand.py(), k0Scand.pz(), massK0s);
        lResoKstar = lResoSecondary + lDecayDaughter_bach;

        // QA plots
        if constexpr (!IsMix) {
          histos.fill(HIST("QA/before/KstarRapidity"), lResoKstar.Rapidity());
          histos.fill(HIST("QA/before/kstarinvmass"), lResoKstar.M());
        }

        if (lResoKstar.Rapidity() > kstarCutCfgs.cKstarMaxRap || lResoKstar.Rapidity() < kstarCutCfgs.cKstarMinRap)
          continue;

        if constexpr (!IsMix) {

          histos.fill(HIST("QA/after/KstarRapidity"), lResoKstar.Rapidity());
          histos.fill(HIST("QA/after/kstarinvmass"), lResoKstar.M());
          histos.fill(HIST("hInvmass_Kstar"), collision.centFT0M(), lResoKstar.Pt(), lResoKstar.M());
          if (helicityCfgs.cBoostKShot) {
            fillInvMass(lResoKstar, collision.centFT0M(), lResoSecondary, lDecayDaughter_bach, IsMix);
          } else {
            fillInvMass(lResoKstar, collision.centFT0M(), lDecayDaughter_bach, lResoSecondary, IsMix);
          }
        } else {

          histos.fill(HIST("hInvmass_KstarME"), collision.centFT0M(), lResoKstar.Pt(), lResoKstar.M());
          fillInvMass(lResoKstar, collision.centFT0M(), lResoSecondary, lDecayDaughter_bach, IsMix);
        }
        if constexpr (!IsMix) {
          if (rotBkgEstCfgs.cFillRotBkg) {
            for (int nrotbkg = 0; nrotbkg < rotBkgEstCfgs.nBkgRotations; nrotbkg++) {
              auto rotangle = o2::constants::math::PI; // If there is only one rotation then it should be pi ):
              if (rotBkgEstCfgs.nBkgRotations > 1) {
                auto anglestart = rotBkgEstCfgs.confMinRot;
                auto angleend = rotBkgEstCfgs.confMaxRot;
                auto anglestep = (angleend - anglestart) / (1.0 * (rotBkgEstCfgs.nBkgRotations - 1));
                rotangle = anglestart + nrotbkg * anglestep;
              }
              histos.fill(HIST("hRotation"), rotangle);
              auto rotpionPx = lDecayDaughter_bach.Px() * std::cos(rotangle) - lDecayDaughter_bach.Py() * std::sin(rotangle);
              auto rotpionPy = lDecayDaughter_bach.Px() * std::sin(rotangle) + lDecayDaughter_bach.Py() * std::cos(rotangle);
              ROOT::Math::PtEtaPhiMVector pionrot;
              pionrot = ROOT::Math::PxPyPzMVector(rotpionPx, rotpionPy, lDecayDaughter_bach.Pz(), massPi);
              chargeKstarrot = pionrot + lResoSecondary;
              if (chargeKstarrot.Rapidity() > kstarCutCfgs.cKstarMaxRap || chargeKstarrot.Rapidity() < kstarCutCfgs.cKstarMinRap)
                continue;
              histos.fill(HIST("hInvmass_KstarRotated"), collision.centFT0M(), chargeKstarrot.Pt(), chargeKstarrot.M());
            }
          }
        }
      } // K0scand
    } // bTrack

    count++;

  } // fillHistograms

  // process data
  void processDataSE(EventCandidates::iterator const& collision,
                     TrackCandidates const& tracks,
                     V0Candidates const& v0s,
                     aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if (!colCuts.isSelected(collision)) // Default event selection
      return;
    histos.fill(HIST("hEvtSelInfo"), 1.5);
    if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 2.5);

    lMultiplicity = getCentrality(collision);
    if (lMultiplicity < eventCutCfgs.cfgEventCentralityMin || lMultiplicity > eventCutCfgs.cfgEventCentralityMax)
      return;
    histos.fill(HIST("hEvtSelInfo"), 3.5);
    if (!collision.isInelGt0())
      return;
    histos.fill(HIST("hEvtSelInfo"), 4.5);
    colCuts.fillQA(collision);
    fillHistograms<false, false>(collision, tracks, v0s);
  }
  PROCESS_SWITCH(Chargedkstaranalysis, processDataSE, "Process Event for data without Partitioning", true);

  SliceCache cache;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningTypeVertexContributor binningOnPositions{{axisCfgs.cfgvtxbins, axisCfgs.cfgmultbins}, true};
  Pair<EventCandidates, TrackCandidates, V0Candidates, BinningTypeVertexContributor> pair{binningOnPositions, nEvtMixing, -1, &cache};

  void processDataME(EventCandidates const& /*collisions*/, TrackCandidates const& /*tracks*/, V0Candidates const& /*V0s*/)
  {
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!colCuts.isSelected(c1)) {
        continue;
      }
      if (!colCuts.isSelected(c2)) {
        continue;
      }

      for (const auto& [t1, t2] : o2::soa::combinations(
             o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        // Here t1 corressponds to bachelor track and t2 corressponds to v0s.
        if (!trackCut(t1))
          continue;
        if (!selectionPIDPion(t1))
          continue;

        auto posDauTrack = t2.template posTrack_as<TrackCandidates>();
        auto negDauTrack = t2.template negTrack_as<TrackCandidates>();
        if (!secondaryCutsCfgs.cfgByPassDauPIDSelection && !selectionPIDPion(posDauTrack)) // Perhaps it's already applied in trackCut (need to check QA plots)
          continue;
        if (!secondaryCutsCfgs.cfgByPassDauPIDSelection && !selectionPIDPion(negDauTrack))
          continue;
        if (!selectionK0s(c2, t2))
          continue;

        ROOT::Math::PxPyPzMVector lResoSecondary, lDecayDaughter_bach, lResoKstar;
        lDecayDaughter_bach = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massPi);
        lResoSecondary = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massK0s);
        lResoKstar = lResoSecondary + lDecayDaughter_bach;

        if (lResoKstar.Rapidity() > kstarCutCfgs.cKstarMaxRap || lResoKstar.Rapidity() < kstarCutCfgs.cKstarMinRap)
          continue;
        histos.fill(HIST("hInvmass_KstarME"), c1.centFT0M(), lResoKstar.Pt(), lResoKstar.M());
      }
    }
  }
  PROCESS_SWITCH(Chargedkstaranalysis, processDataME, "Process Event for data without Partitioning", true);

  // process MC reconstructed level
  void processMC(MCEventCandidates::iterator const& collision,
                 MCTrackCandidates const& tracks,
                 MCV0Candidates const& v0s)
  {

    //    histos.fill(HIST("QAMC/hEvent"), 1.0);

    fillHistograms<true, false>(collision, tracks, v0s);
  }
  PROCESS_SWITCH(Chargedkstaranalysis, processMC, "Process Event for MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Chargedkstaranalysis>(cfgc)};
}
