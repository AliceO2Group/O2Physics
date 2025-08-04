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

/// \file chk892_pp.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>

#include <TH1F.h>
#include <TH1D.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TDatabasePDG.h> // FIXME
#include <TPDGCode.h>     // FIXME

#include <vector>
#include <cmath>
#include <array>
#include <cstdlib>
#include <chrono>
#include <string>

#include "TRandom3.h"
#include "TF1.h"
#include "TVector2.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"
#include "DCAFitter/DCAFitterN.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/RecoDecay.h"

#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/collisionCuts.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct chk892_pp {
  enum binType : unsigned int {
    kKstarP = 0,
    kKstarN,
    kKstarP_Mix,
    kKstarN_Mix,
    kKstarP_GenINEL10,
    kKstarN_GenINEL10,
    kKstarP_GenINELgt10,
    kKstarN_GenINELgt10,
    kKstarP_GenTrig10,
    kKstarN_GenTrig10,
    kKstarP_GenEvtSel,
    kKstarN_GenEvtSel,
    kKstarP_Rec,
    kKstarN_Rec,
    kTYEnd
  };

  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  //  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::Mults>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using V0Candidates = aod::V0Datas;

  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Join<TrackCandidates, aod::McTrackLabels>;
  using MCV0Candidates = soa::Join<V0Candidates, aod::McV0Labels>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  // Configurables
  ConfigurableAxis cfgBinsPt{"cfgBinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis cfgBinsPtQA{"cfgBinsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis cfgBinsCent{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  ConfigurableAxis cfgBinsVtxZ{"cfgBinsVtxZ", {VARIABLE_WIDTH, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "Binning of the z-vertex axis"};
  Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<int> ConfEvtOccupancyInTimeRangeMax{"ConfEvtOccupancyInTimeRangeMax", -1, "Evt sel: maximum track occupancy"};
  Configurable<int> ConfEvtOccupancyInTimeRangeMin{"ConfEvtOccupancyInTimeRangeMin", -1, "Evt sel: minimum track occupancy"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", false, "Evt sel: check for trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", true, "Evt sel: check for offline selection"};
  Configurable<bool> ConfEvtTriggerTVXSel{"ConfEvtTriggerTVXSel", false, "Evt sel: triggerTVX selection (MB)"};
  Configurable<bool> ConfEvtTFBorderCut{"ConfEvtTFBorderCut", false, "Evt sel: apply TF border cut"};
  Configurable<bool> ConfEvtUseITSTPCvertex{"ConfEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
  Configurable<bool> ConfEvtZvertexTimedifference{"ConfEvtZvertexTimedifference", true, "Evt sel: apply Z-vertex time difference"};
  Configurable<bool> ConfEvtPileupRejection{"ConfEvtPileupRejection", true, "Evt sel: apply pileup rejection"};
  Configurable<bool> ConfEvtNoITSROBorderCut{"ConfEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
  Configurable<bool> ConfincludeCentralityMC{"ConfincludeCentralityMC", false, "Include centrality in MC"};
  Configurable<bool> ConfEvtCollInTimeRangeStandard{"ConfEvtCollInTimeRangeStandard", true, "Evt sel: apply NoCollInTimeRangeStandard"};

  /// Track selections
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  Configurable<float> cMaxEtacut{"cMaxEtacut", 0.8, "Track maximum eta cut"};

  /*
  // Cuts from polarization analysis
  Configurable<bool> cfgOccupancySel{"cfgOccupancySel", false, "Occupancy selection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", -100, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<bool> cfgNCollinTR{"cfgNCollinTR", false, "Additional selection for the number of coll in time range"};
*/
  Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  // DCAr to PV
  Configurable<float> cMaxbDCArToPVcut{"cMaxbDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<float> cMaxbDCAzToPVcut{"cMaxbDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};

  /// PID Selections, pion
  Configurable<bool> cTPConly{"cTPConly", true, "Use only TPC for PID"};                                    // bool
  Configurable<float> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};               // TPC
  Configurable<float> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};               // TOF
  Configurable<float> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"}; // Combined
  Configurable<bool> cTOFVeto{"cTOFVeto", true, "TOF Veto, if false, TOF is nessessary for PID selection"}; // TOF Veto

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

  // Secondary Selection
  Configurable<bool> cfgReturnFlag{"boolReturnFlag", false, "Return Flag for debugging"};
  Configurable<bool> cSecondaryRequire{"bool", true, "Secondary cuts on/off"};
  Configurable<bool> cSecondaryArmenterosCut{"boolArmenterosCut", true, "cut on Armenteros-Podolanski graph"};

  Configurable<bool> cfgByPassDauPIDSelection{"cfgByPassDauPIDSelection", true, "Bypass Daughters PID selection"};
  Configurable<float> cSecondaryDauDCAMax{"cSecondaryDauDCAMax", 1., "Maximum DCA Secondary daughters to PV"};
  Configurable<float> cSecondaryDauPosDCAtoPVMin{"cSecondaryDauPosDCAtoPVMin", 0.0, "Minimum DCA Secondary positive daughters to PV"};
  Configurable<float> cSecondaryDauNegDCAtoPVMin{"cSecondaryDauNegDCAtoPVMin", 0.0, "Minimum DCA Secondary negative daughters to PV"};

  Configurable<float> cSecondaryPtMin{"cSecondaryPtMin", 0.f, "Minimum transverse momentum of Secondary"};
  Configurable<float> cSecondaryRapidityMax{"cSecondaryRapidityMax", 0.5, "Maximum rapidity of Secondary"};
  Configurable<float> cSecondaryRadiusMin{"cSecondaryRadiusMin", 1.2, "Minimum transverse radius of Secondary"};
  Configurable<float> cSecondaryCosPAMin{"cSecondaryCosPAMin", 0.995, "Mininum cosine pointing angle of Secondary"};
  Configurable<float> cSecondaryDCAtoPVMax{"cSecondaryDCAtoPVMax", 0.3, "Maximum DCA Secondary to PV"};
  Configurable<float> cSecondaryProperLifetimeMax{"cSecondaryProperLifetimeMax", 20, "Maximum Secondary Lifetime"};
  Configurable<float> cSecondaryparamArmenterosCut{"paramArmenterosCut", 0.2, "parameter for Armenteros Cut"};
  Configurable<float> cSecondaryMassWindow{"cSecondaryMassWindow", 0.075, "Secondary inv mass selciton window"};

  // K* selection
  Configurable<float> cKstarMaxRap{"cKstarMaxRap", 0.5, "Kstar maximum rapidity"};
  Configurable<float> cKstarMinRap{"cKstarMinRap", -0.5, "Kstar minimum rapidity"};

  float centrality;

  // PDG code
  int kPDGK0s = 310;
  int kPDGK0 = 311;
  int kKstarPlus = 323;
  int kPiPlus = 211;

  void init(o2::framework::InitContext&)
  {
    centrality = -999;

    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtOfflineCheck, /*checkRun3*/ true, /*triggerTVXsel*/ false, ConfEvtOccupancyInTimeRangeMax, ConfEvtOccupancyInTimeRangeMin);
    colCuts.init(&histos);
    colCuts.setTriggerTVX(ConfEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(ConfEvtTFBorderCut);
    colCuts.setApplyITSTPCvertex(ConfEvtUseITSTPCvertex);
    colCuts.setApplyZvertexTimedifference(ConfEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(ConfEvtPileupRejection);
    colCuts.setApplyNoITSROBorderCut(ConfEvtNoITSROBorderCut);
    colCuts.setApplyCollInTimeRangeStandard(ConfEvtCollInTimeRangeStandard);

    AxisSpec centAxis = {cfgBinsCent, "T0M (%)"};
    AxisSpec vtxzAxis = {cfgBinsVtxZ, "Z Vertex (cm)"};
    AxisSpec epAxis = {100, -1.0 * constants::math::PI, constants::math::PI};
    AxisSpec epresAxis = {100, -1.02, 1.02};
    AxisSpec ptAxis = {cfgBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {cfgBinsPtQA, "#it{p}_{T} (GeV/#it{c})"};
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

    // THnSparse
    AxisSpec axisType = {binType::kTYEnd, 0, binType::kTYEnd, "Type of bin with charge and mix"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};

    histos.add("QA/K0sCutCheck", "Check K0s cut", HistType::kTH1D, {AxisSpec{12, -0.5, 11.5, "Check"}});

    histos.add("QA/before/CentDist", "Centrality distribution", {HistType::kTH1D, {centAxis}});
    histos.add("QA/before/VtxZ", "Centrality distribution", {HistType::kTH1D, {vtxzAxis}});
    histos.add("QA/before/hEvent", "Number of Events", HistType::kTH1F, {{1, 0.5, 1.5}});

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

    // Mass QA (quick check)
    histos.add("QA/before/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
    histos.add("QA/before/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

    histos.add("QA/after/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});
    histos.add("QA/after/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
    histos.add("QA/after/kstarinvmass_Mix", "Invariant mass of unlike-sign chK(892) from mixed event", HistType::kTH1D, {invMassAxisReso});

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
      histos.add("QAMC/KstarPairAsym", "Pair asymmetry of chK(892)", HistType::kTH1D, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
      histos.add("QAMC/KstarRapidity", "Rapidity distribution of chK(892)", HistType::kTH1D, {yAxis});

      histos.add("QAMC/kstarinvmass", "Invariant mass of unlike-sign chK(892)", HistType::kTH1D, {invMassAxisReso});
      histos.add("QAMC/kstarinvmass_noKstar", "Invariant mass of unlike-sign no chK(892)", HistType::kTH1D, {invMassAxisReso});

      histos.add("hInvmass_Kstar_MC", "Invariant mass of unlike chK(892)", HistType::kTHnSparseD, {axisType, centAxis, ptAxis, invMassAxisReso});

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

  template <typename CollisionType>
  float GetCentrality(CollisionType const& collision)
  {
    if (cfgCentEst == 1) {
      return collision.multFT0C();
    } else if (cfgCentEst == 2) {
      return collision.multFT0M();
    } else {
      return -999;
    }
  }

  template <typename DetNameType>
  int GetDetId(DetNameType const& name)
  {
    LOGF(info, "GetDetID running");
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
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.eta()) > cMaxEtacut)
      return false;
    if (track.itsNCls() < cfgITScluster)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < cfgRatioTPCRowsOverFindableCls)
      return false;
    if (track.itsChi2NCl() >= cfgITSChi2NCl)
      return false;
    if (track.tpcChi2NCl() >= cfgTPCChi2NCl)
      return false;
    if (cfgHasITS && !track.hasITS())
      return false;
    if (cfgHasTPC && !track.hasTPC())
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (std::abs(track.dcaXY()) > cMaxbDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxbDCAzToPVcut)
      return false;
    return true;
  }

  // PID selection tools
  template <typename TrackType>
  bool selectionPIDPion(TrackType const& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (cTPConly) {

      if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
        tpcPIDPassed = true;
      } else {
        return false;
      }
      tofPIDPassed = true;

    } else {

      if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
        tpcPIDPassed = true;
      } else {
        return false;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion) {
          tofPIDPassed = true;
        }
        if ((nsigmaCutCombinedPion > 0) && (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
          tofPIDPassed = true;
        }
      } else {
        if (!cTOFVeto) {
          return false;
        }
        tofPIDPassed = true;
      }
    }

    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <typename CollisionType, typename K0sType>
  bool selectionK0s(CollisionType const& collision, K0sType const& candidate)
  {
    auto DauDCA = candidate.dcaV0daughters();
    auto DauPosDCAtoPV = candidate.dcapostopv();
    auto DauNegDCAtoPV = candidate.dcanegtopv();
    auto pT = candidate.pt();
    auto Rapidity = candidate.yK0Short();
    auto Radius = candidate.v0radius();
    auto DCAtoPV = candidate.dcav0topv();
    auto CPA = candidate.v0cosPA();
    auto PropTauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassK0Short;
    auto mK0s = candidate.mK0Short();

    if (cfgReturnFlag) {
      bool returnFlag = true;

      if (cSecondaryRequire) {
        histos.fill(HIST("QA/K0sCutCheck"), 0);
        if (DauDCA > cSecondaryDauDCAMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 1);
          returnFlag = false;
        }
        if (DauPosDCAtoPV < cSecondaryDauPosDCAtoPVMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 2);
          returnFlag = false;
        }
        if (DauNegDCAtoPV < cSecondaryDauNegDCAtoPVMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 3);
          returnFlag = false;
        }
        if (pT < cSecondaryPtMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 4);
          returnFlag = false;
        }
        if (Rapidity > cSecondaryRapidityMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 5);
          returnFlag = false;
        }
        if (Radius < cSecondaryRadiusMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 6);
          returnFlag = false;
        }
        if (DCAtoPV > cSecondaryDCAtoPVMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 7);
          returnFlag = false;
        }
        if (CPA < cSecondaryCosPAMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 8);
          returnFlag = false;
        }
        if (PropTauK0s > cSecondaryProperLifetimeMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 9);
          returnFlag = false;
        }
        if (candidate.qtarm() < cSecondaryparamArmenterosCut * TMath::Abs(candidate.alpha())) {
          histos.fill(HIST("QA/K0sCutCheck"), 11);
          returnFlag = false;
        }
        if (fabs(mK0s - MassK0Short) > cSecondaryMassWindow) {
          histos.fill(HIST("QA/K0sCutCheck"), 10);
          returnFlag = false;
        }

        return returnFlag;

      } else {
        if (fabs(mK0s - MassK0Short) > cSecondaryMassWindow) {
          histos.fill(HIST("QA/K0sCutCheck"), 10);
          returnFlag = false;
        }

        return returnFlag;
      }

    } else {
      if (cSecondaryRequire) {

        histos.fill(HIST("QA/K0sCutCheck"), 0);
        if (DauDCA > cSecondaryDauDCAMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 1);
          return false;
        }
        if (DauPosDCAtoPV < cSecondaryDauPosDCAtoPVMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 2);
          return false;
        }
        if (DauNegDCAtoPV < cSecondaryDauNegDCAtoPVMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 3);
          return false;
        }
        if (pT < cSecondaryPtMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 4);
          return false;
        }
        if (Rapidity > cSecondaryRapidityMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 5);
          return false;
        }
        if (Radius < cSecondaryRadiusMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 6);
          return false;
        }
        if (DCAtoPV > cSecondaryDCAtoPVMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 7);
          return false;
        }
        if (CPA < cSecondaryCosPAMin) {
          histos.fill(HIST("QA/K0sCutCheck"), 8);
          return false;
        }
        if (PropTauK0s > cSecondaryProperLifetimeMax) {
          histos.fill(HIST("QA/K0sCutCheck"), 9);
          return false;
        }
        if (candidate.qtarm() < cSecondaryparamArmenterosCut * TMath::Abs(candidate.alpha())) {
          histos.fill(HIST("QA/K0sCutCheck"), 11);
          return false;
        }
        if (fabs(mK0s - MassK0Short) > cSecondaryMassWindow) {
          histos.fill(HIST("QA/K0sCutCheck"), 10);
          return false;
        }
        return true;

      } else {
        if (fabs(mK0s - MassK0Short) > cSecondaryMassWindow) {
          histos.fill(HIST("QA/K0sCutCheck"), 10);
          return false;
        }
        return true;
      }
    }
  } // selectionK0s

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }

  template <typename TrackTemplate, typename V0Template>
  bool isTrueKstar(const TrackTemplate& bTrack, const V0Template& K0scand)
  {
    if (abs(bTrack.PDGCode()) != kPiPlus) // Are you pion?
      return false;
    if (abs(K0scand.PDGCode()) != kPDGK0s) // Are you K0s?
      return false;

    auto motherbTrack = bTrack.template mothers_as<aod::McParticles>();
    auto motherkV0 = K0scand.template mothers_as<aod::McParticles>();

    // Check bTrack first
    if (abs(motherbTrack.pdgCode()) != kKstarPlus) // Are you charged Kstar's daughter?
      return false;                                // Apply first since it's more restrictive

    if (abs(motherkV0.pdgCode()) != 310) // Is it K0s?
      return false;
    // Check if K0s's mother is K0 (311)
    auto motherK0 = motherkV0.template mothers_as<aod::McParticles>();
    if (abs(motherK0.pdgCode()) != 311)
      return false;

    // Check if K0's mother is Kstar (323)
    auto motherKstar = motherK0.template mothers_as<aod::McParticles>();
    if (abs(motherKstar.pdgCode()) != 323)
      return false;

    // Check if bTrack and K0 have the same mother (global index)
    if (motherbTrack.globalIndex() != motherK0.globalIndex())
      return false;

    return true;
  }

  int count = 0;

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType, typename TracksTypeK0s>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksTypeK0s& dTracks2)
  {
    histos.fill(HIST("QA/before/CentDist"), centrality);

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResoSecondary, lDecayDaughter_bach, lResoKstar;
    std::vector<int> trackIndicies = {};
    std::vector<int> k0sIndicies = {};

    for (auto& bTrack : dTracks1) {
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

    for (auto& K0scand : dTracks2) {
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
      auto trkkPropTau = K0scand.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassK0Short;
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
        histos.fill(HIST("QA/before/hPropTauSecondary"), trkkPropTau);
        histos.fill(HIST("QA/before/hInvmassSecondary"), trkkMass);
      }

      // if (!trackCut(posDauTrack) || !trackCut(negDauTrack)) // Too tight cut for K0s daugthers
      //   continue;
      if (!cfgByPassDauPIDSelection && !selectionPIDPion(posDauTrack)) // Perhaps it's already applied in trackCut (need to check QA plots)
        continue;
      if (!cfgByPassDauPIDSelection && !selectionPIDPion(negDauTrack))
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
        histos.fill(HIST("QA/after/hPropTauSecondary"), trkkPropTau);
        histos.fill(HIST("QA/after/hInvmassSecondary"), trkkMass);
      }
      k0sIndicies.push_back(K0scand.index());
    }

    for (auto& trackIndex : trackIndicies) {
      for (auto& k0sIndex : k0sIndicies) {
        auto bTrack = dTracks1.rawIteratorAt(trackIndex);
        auto K0scand = dTracks2.rawIteratorAt(k0sIndex);

        lDecayDaughter_bach.SetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), MassPionCharged);
        lResoSecondary.SetXYZM(K0scand.px(), K0scand.py(), K0scand.pz(), MassK0Short);
        lResoKstar = lResoSecondary + lDecayDaughter_bach;

        // QA plots
        if constexpr (!IsMix) {
          histos.fill(HIST("QA/before/KstarRapidity"), lResoKstar.Rapidity());
          histos.fill(HIST("QA/before/kstarinvmass"), lResoKstar.M());
        }

        if (lResoKstar.Rapidity() > cKstarMaxRap || lResoKstar.Rapidity() < cKstarMinRap)
          continue;

        if constexpr (!IsMix) {
          unsigned int typeKstar = bTrack.sign() > 0 ? binType::kKstarP : binType::kKstarN;

          histos.fill(HIST("QA/after/KstarRapidity"), lResoKstar.Rapidity());
          histos.fill(HIST("QA/after/kstarinvmass"), lResoKstar.M());
          histos.fill(HIST("hInvmass_Kstar"), typeKstar, centrality, lResoKstar.Pt(), lResoKstar.M());

        } // IsMix
      } // K0scand
    } // bTrack

    count++;

  } // fillHistograms

  // process data
  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks,
                   V0Candidates const& v0s,
                   aod::BCsWithTimestamps const&)
  {
    if (!colCuts.isSelected(collision)) // Default event selection
      return;
    colCuts.fillQA(collision);
    centrality = GetCentrality(collision);

    fillHistograms<false, false>(collision, tracks, v0s); // second order
  }
  PROCESS_SWITCH(chk892_pp, processData, "Process Event for data without Partitioning", true);

  // process MC reconstructed level
  void processMC(EventCandidates::iterator const& collision,
                 MCTrackCandidates const& tracks,
                 MCV0Candidates const& v0s)
  {

    histos.fill(HIST("QAMC/hEvent"), 1.0);

    fillHistograms<true, false>(collision, tracks, v0s);
  }
  PROCESS_SWITCH(chk892_pp, processMC, "Process Event for MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<chk892_pp>(cfgc, TaskName{"lf-chk892pp"})};
}
