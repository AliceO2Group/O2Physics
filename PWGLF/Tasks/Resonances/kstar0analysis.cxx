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

/// \file kstar0analysis.cxx
/// \brief the reconstruction of k*0(892) resonance analysis using K\pi decay channel
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

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
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TRandom.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::aod::rctsel;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using LorentzVectorPtEtaPhiMass = ROOT::Math::PtEtaPhiMVector;

enum TrackSelectionType {
  AllTracks = 0,
  GlobalTracks,
  GlobalTracksWoPtEta,
  GlobalTracksWoDCA,
  QualityTracks,
  InAcceptanceTracks,
};

enum PIDCutType {
  SquareType = 1,
  CircularType,
};

struct Kstar0analysis {
  // Define slice per collision
  Preslice<Tracks> perCollision = o2::aod::track::collisionId;
  SliceCache cache;
  Preslice<McParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;
  SliceCache cacheMC;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<framework::O2DatabasePDG> pdg;
  RCTFlagsChecker rctChecker;

  struct : ConfigurableGroup {
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.0f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", true, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtNoTFBorderCut{"cfgEvtNoTFBorderCut", true, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtNoITSROFrameBorderCut{"cfgEvtNoITSROFrameBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtIsRCTFlagpassed{"cfgEvtIsRCTFlagpassed", false, "Evt sel: apply RCT flag selection"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
    Configurable<bool> cfgEvtSel8{"cfgEvtSel8", false, "Evt Sel 8 check for offline selection"};
    Configurable<bool> cfgEvtIsINELgt0{"cfgEvtIsINELgt0", false, "Evt sel: apply INEL>0 selection"};
  } configEvents;

  struct : ConfigurableGroup {
    // Pre-selection Track cuts
    Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
    Configurable<float> cMinPtcut{"cMinPtcut", 0.15f, "Minimal pT for tracks"};
    Configurable<float> cMinTPCNClsFound{"cMinTPCNClsFound", 120, "minimum TPCNClsFound value for good track"};
    Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
    Configurable<float> cfgCutRapidity{"cfgCutRapidity", 0.5f, "rapidity range for particles"};
    Configurable<int> cfgMinCrossedRows{"cfgMinCrossedRows", 70, "min crossed rows for good track"};

    // DCA Selections
    // DCAr to PV
    Configurable<float> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1f, "Track DCAr cut to PV Maximum"};
    // DCAz to PV
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1f, "Track DCAz cut to PV Maximum"};

    // Track selections
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
    Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
    Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
    Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
    Configurable<bool> cTPCNClsFound{"cTPCNClsFound", false, "Switch to turn on/off TPCNClsFound cut"};
    Configurable<bool> cDCAr7SigCut{"cDCAr7SigCut", false, "Track DCAr 7 Sigma cut to PV Maximum"};
  } configTracks;

  struct : ConfigurableGroup {
    /// PID Selections
    Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 4.0f, "pidnSigma Cut for pre-selection of tracks"};
    Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};                       // By pass TOF PID selection
    Configurable<int> cPIDcutType{"cPIDcutType", 2, "cPIDcutType = 1 for square cut, 2 for circular cut"}; // By pass TOF PID selection

    // Kaon
    Configurable<std::vector<float>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {0.5f}, "pT intervals for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.0f}, "pT intervals for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> kaonTPCTOFCombinedpTintv{"kaonTPCTOFCombinedpTintv", {999.0f}, "pT intervals for Kaon TPC-TOF PID cuts"};
    Configurable<std::vector<float>> kaonTPCTOFCombinedPIDcuts{"kaonTPCTOFCombinedPIDcuts", {2}, "nSigma list for Kaon TPC-TOF PID cuts"};

    // Pion
    Configurable<std::vector<float>> pionTPCPIDpTintv{"pionTPCPIDpTintv", {0.9f}, "pT intervals for Pion TPC PID cuts"};
    Configurable<std::vector<float>> pionTPCPIDcuts{"pionTPCPIDcuts", {2}, "nSigma list for Pion TPC PID cuts"};
    Configurable<std::vector<float>> pionTOFPIDpTintv{"pionTOFPIDpTintv", {999.0f}, "pT intervals for Pion TOF PID cuts"};
    Configurable<std::vector<float>> pionTOFPIDcuts{"pionTOFPIDcuts", {2}, "nSigma list for Pion TOF PID cuts"};
    Configurable<std::vector<float>> pionTPCTOFCombinedpTintv{"pionTPCTOFCombinedpTintv", {999.0f}, "pT intervals for Pion TPC-TOF PID cuts"};
    Configurable<std::vector<float>> pionTPCTOFCombinedPIDcuts{"pionTPCTOFCombinedPIDcuts", {2}, "nSigma list for Pion TPC-TOF PID cuts"};
  } configPID;

  struct : ConfigurableGroup {
    /// Event Mixing
    Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgMultPercentileBins{"cfgMultPercentileBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - multiplicity"};

    // Rotational background
    Configurable<int> rotationalcut{"rotationalcut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
    Configurable<int> cNofRotations{"cNofRotations", 3, "Number of random rotations in the rotational background"};
    Configurable<bool> cfgRotPi{"cfgRotPi", true, "rotate Pion"};
  } configBkg;

  // Additional purity check
  Configurable<bool> crejectProton{"crejectProton", false, "Switch to turn on/off proton contamination"};
  Configurable<bool> cUseOpeningAngleCut{"cUseOpeningAngleCut", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<float> cMinOpeningAngle{"cMinOpeningAngle", 1.4f, "Minimum opening angle between daughters"};
  Configurable<float> cMaxOpeningAngle{"cMaxOpeningAngle", 2.4f, "Maximum opening angle between daughters"};
  Configurable<bool> cfgUseDeltaEtaPhiCuts{"cfgUseDeltaEtaPhiCuts", false, "Switch to turn on/off delta eta and delta phi cuts"};
  Configurable<bool> cfgUseDaughterEtaCutMC{"cfgUseDaughterEtaCutMC", false, "Switch to turn on/off eta cuts for daughters in MC"};

  // MC selection cut
  Configurable<bool> cUseRapcutMC{"cUseRapcutMC", true, "Use rapidity cut for MC"};

  // cuts on mother
  Configurable<bool> cfgUseCutsOnMother{"cfgUseCutsOnMother", false, "Enable additional cuts on mother"};
  Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 10.0f, "Maximum pt of mother cut"};
  Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 3.0f, "Maximum Minv of mother cut"};
  Configurable<float> cMaxDeltaEtaCut{"cMaxDeltaEtaCut", 0.7f, "Maximum deltaEta between daughters"};
  Configurable<float> cMaxDeltaPhiCut{"cMaxDeltaPhiCut", 1.5f, "Maximum deltaPhi between daughters"};

  // switches
  Configurable<bool> cFillMultQA{"cFillMultQA", false, "Turn on/off additional QA plots"};
  Configurable<bool> cFillTrackQA{"cFillTrackQA", false, "Turn on/off additional QA plots"};
  Configurable<bool> cFilladditionalQAeventPlots{"cFilladditionalQAeventPlots", false, "Additional QA event plots"};
  Configurable<bool> cFilladditionalMEPlots{"cFilladditionalMEPlots", false, "Additional Mixed event plots"};
  Configurable<bool> cFilldeltaEtaPhiPlots{"cFilldeltaEtaPhiPlots", false, "Enable additional cuts on daughters"};
  Configurable<bool> cFill1DQAs{"cFill1DQAs", false, "Invariant mass 1D"};
  Configurable<int> centEstimator{"centEstimator", 0, "Select centrality estimator: 0 - FT0M, 1 - FT0A, 2 - FT0C"};

  TRandom* rn = new TRandom();

  // Pre-filters for efficient process
  Filter zVtxFilter = (nabs(o2::aod::collision::posZ) <= configEvents.cfgEvtZvtx);
  Filter collisionFilterMC = nabs(aod::mccollision::posZ) <= configEvents.cfgEvtZvtx;

  Filter acceptanceFilter = (nabs(aod::track::eta) < configTracks.cfgCutEta && nabs(aod::track::pt) > configTracks.cMinPtcut) &&
                            (nabs(aod::track::dcaXY) < configTracks.cMaxDCArToPVcut) && (nabs(aod::track::dcaZ) < configTracks.cMaxDCAzToPVcut);

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::FullTracks, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>>;

  // for MC reco
  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Filtered<soa::Join<TrackCandidates, aod::McTrackLabels>>;

  /// Figures
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.12f, 0.14f, 0.16f, 0.18f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.45f, 0.5f, 0.55f, 0.6f, 0.65f, 0.7f, 0.75f, 0.8f, 0.85f, 0.9f, 0.95f, 1.0f, 1.1f, 1.2f, 1.25f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.75f, 1.8f, 1.9f, 2.0f, 2.1f, 2.2f, 2.3f, 2.4f, 2.5f, 2.6f, 2.7f, 2.8f, 2.9f, 3.0f, 3.1f, 3.2f, 3.3f, 3.4f, 3.6f, 3.7f, 3.8f, 3.9f, 4.0f, 4.1f, 4.2f, 4.5f, 4.6f, 4.8f, 4.9f, 5.0f, 5.5f, 5.6f, 6.0f, 6.4f, 6.5f, 7.0f, 7.2f, 8.0f, 9.0f, 9.5f, 9.6f, 10.0f, 11.0f, 11.5f, 12.0f, 13.0f, 14.0f, 14.4f, 15.0f, 16.0f, 18.0f, 19.2f, 20.0f}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.2f, 4.4f, 4.6f, 4.8f, 5.0f, 5.2f, 5.4f, 5.6f, 5.8f, 6.0f, 6.2f, 6.4f, 6.6f, 6.8f, 7.0f, 7.2f, 7.4f, 7.6f, 7.8f, 8.0f, 8.2f, 8.4f, 8.6f, 8.8f, 9.0f, 9.2f, 9.4f, 9.6f, 9.8f, 10.0f, 10.2f, 10.4f, 10.6f, 10.8f, 11.0f, 11.2f, 11.4f, 11.6f, 11.8f, 12.0f, 12.2f, 12.4f, 12.6f, 12.8f, 13.0f, 13.2f, 13.4f, 13.6f, 13.8f, 14.0f, 14.2f, 14.4f, 14.6f, 14.8f, 15.0f, 15.2f, 15.4f, 15.6f, 15.8f, 16.0f, 16.2f, 16.4f, 16.6f, 16.8f, 17.0f, 17.2f, 17.4f, 17.6f, 17.8f, 18.0f, 18.2f, 18.4f, 18.6f, 18.8f, 19.0f, 19.2f, 19.4f, 19.6f, 19.8f, 20.0f}, "Binning of the pT axis"};
  ConfigurableAxis binsEta{"binsEta", {150, -1.5f, 1.5f}, ""};
  ConfigurableAxis binsMass{"binsMass", {300, 0.6f, 1.2f}, "Invariant Mass (GeV/#it{c}^2)"};
  ConfigurableAxis binsMult{"binsMult", {105, 0.0f, 105.0f}, "mult_{FT0M}"};
  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2f, 0.2f}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2f, 0.2f}, ""};
  ConfigurableAxis binsTPCXrows{"binsTPCXrows", {100, 60, 160}, ""};
  ConfigurableAxis binsnSigma{"binsnSigma", {130, -6.5f, 6.5f}, ""};
  ConfigurableAxis binsnTPCSignal{"binsnTPCSignal", {1000, 0, 1000}, ""};
  ConfigurableAxis binsEtaPhi{"binsEtaPhi", {350, -3.5f, 3.5f}, ""};

  void init(o2::framework::InitContext&)
  {
    rctChecker.init(configEvents.cfgEvtRCTFlagCheckerLabel, configEvents.cfgEvtRCTFlagCheckerZDCCheck, configEvents.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    // axes
    AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtQA{binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{binsEta, "#eta"};
    AxisSpec axisRap{binsEta, "#it{y}"};
    AxisSpec axisMassk892{binsMass, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec axisMult{binsMult, "mult_{FT0M}"};
    AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};
    AxisSpec axisTPCXrow{binsTPCXrows, "#Xrows_{TPC}"};
    AxisSpec axisPIDQA{binsnSigma, "#sigma"};
    AxisSpec axisTPCSignal{binsnTPCSignal, ""};
    AxisSpec axisEtaPhi{binsEtaPhi, ""};
    AxisSpec axisPhi{350, 0, 7, "#Phi"};
    AxisSpec axisMultMix{configBkg.cfgMultPercentileBins, "Multiplicity Percentile"};
    AxisSpec axisVtxMix{configBkg.cfgVtxBins, "Vertex Z (cm)"};

    histos.add("CollCutCounts", "No. of event after cuts", kTH1I, {{10, 0, 10}});
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(1, "All Events");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(2, "|Vz| < cut");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(5, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(6, "rctChecker");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(7, "sel8");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(8, "IsINELgt0");
    histos.get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(9, "All Passed Events");

    histos.add("Event/posZ", "; vtx_{z} (cm); Entries", HistType::kTH1F, {{250, -12.5, 12.5}});
    histos.add("Event/centFT0M", "; FT0M Percentile; Entries", HistType::kTH1F, {{110, 0, 110}});

    if (cFilladditionalQAeventPlots) {
      // event histograms
      if (doprocessData) {
        histos.add("QAevent/hPairsCounterSameE", "total valid no. of pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
        histos.add("QAevent/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0, 1000.0}});
      }
      // Gen on Mixed event
      if (doprocessME) {

        // Histograms for Mixed Event Pool characteristics
        histos.add("QAevent/hMixPool_VtxZ", "Mixed Event Pool: Vertex Z;Vtx Z (cm);Counts", HistType::kTH1F, {axisVtxMix});
        histos.add("QAevent/hMixPool_Multiplicity", "Mixed Event Pool: Multiplicity;Multiplicity;Counts", HistType::kTH1F, {axisMultMix});
        histos.add("QAevent/hMixPool_VtxZ_vs_Multiplicity", "Mixed Event Pool: Vertex Z vs Multiplicity;Counts", HistType::kTH2F, {axisVtxMix, axisMultMix});

        histos.add("QAevent/hPairsCounterMixedE", "total valid no. of pairs mixedE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
        histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15.0f, 15.0f}});
        histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
        histos.add("QAevent/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      }
    }

    if (doprocessData) {
      // Track QA before cuts
      //  --- Track
      if (cFillTrackQA) {
        histos.add("QA/QAbefore/Track/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for all Kaon;{#sigma_{TOF}^{Kaon}};{#sigma_{TPC}^{Kaon}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TOF_Nsigma_ka_all", "TOF NSigma for all Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPC_Nsigma_ka_all", "TPC NSigma for all Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPConly_Nsigma_ka_all", "TPC NSigma for all Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTH2F, {axisPt, axisPIDQA}});

        histos.add("QA/QAbefore/Track/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for all Pion;{#sigma_{TOF}^{Pion}};{#sigma_{TPC}^{Pion}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TOF_Nsigma_pi_all", "TOF NSigma for all Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Pion}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPC_Nsigma_pi_all", "TPC NSigma for all Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Pion}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPConly_Nsigma_pi_all", "TPC NSigma for all Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Pion}};", {HistType::kTH2F, {axisPt, axisPIDQA}});

        histos.add("QA/QAbefore/Track/dcaZ_all", "DCA_{Z} distribution of all Tracks; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
        histos.add("QA/QAbefore/Track/dcaXY_all", "DCA_{XY} momentum distribution of all Tracks; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
        histos.add("QA/QAbefore/Track/TPC_CR_all", "# TPC Xrows distribution of all Tracks; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
        histos.add("QA/QAbefore/Track/pT_all", "pT distribution of all Tracks; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
        histos.add("QA/QAbefore/Track/eta_all", "#eta distribution of all Tracks; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      }
      if (cFillMultQA) {
        // Multiplicity correlation calibrations
        histos.add("MultCalib/centGloPVpi", "Centrality vs Global-Tracks", kTHnSparseF, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}, {500, 0, 5000, "PV tracks"}});
        histos.add("MultCalib/centGloPVka", "Centrality vs Global-Tracks", kTHnSparseF, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}, {500, 0, 5000, "PV tracks"}});
      }

      // PID QA after cuts
      //  --- Kaon
      histos.add("QA/QAafter/Kaon/TOF_TPC_Map_ka_selected", "TOF + TPC Combined PID for selected Kaon;{#sigma_{TOF}^{Kaon}};{#sigma_{TPC}^{Kaon}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/TOF_Nsigma_ka_selected", "TOF NSigma for selected Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_selected", "TPC NSigma for selected Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly_selected", "TPC NSigma for selected Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTH2F, {axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/dcaZ_selected", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Kaon/dcaXY_selected", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Kaon/TPC_CR_selected", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Kaon/pT_selected", "pT distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Kaon/eta_selected", "#eta distribution of selected Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Kaon/TPC_Signal_ka_selected", "TPC Signal for selected Kaon;#it{p} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});
      histos.add("QA/QAafter/Kaon/TPCnclusterPhika_selected", "TPC ncluster vs phi for selected Kaon", kTHnSparseF, {{160, 0, 160, "TPC nCluster"}, {63, 0.0f, 6.28f, "#phi"}});

      //  --- Pion
      histos.add("QA/QAafter/Pion/TOF_TPC_Map_pi_selected", "TOF + TPC Combined PID for selected Pion;{#sigma_{TOF}^{Pion}};{#sigma_{TPC}^{Pion}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
      histos.add("QA/QAafter/Pion/TOF_Nsigma_pi_selected", "TOF NSigma for selected Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Pion}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Pion/TPC_Nsigma_pi_selected", "TPC NSigma for selected Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Pion}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Pion/TPC_Nsigma_pi_TPConly_selected", "TPC NSigma for selected Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Pion}};", {HistType::kTH2F, {axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Pion/dcaZ_selected", "DCA_{Z} distribution of selected Pions; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm);", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Pion/dcaXY_selected", "DCA_{XY} momentum distribution of selected Pions; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Pion/TPC_CR_selected", "# TPC Xrows distribution of selected Pions; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Pion/pT_selected", "pT distribution of selected Pions; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Pion/eta_selected", "#eta distribution of selected Pions; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Pion/TPC_Signal_pi_selected", "TPC Signal for selected Pion;#it{p} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});
      histos.add("QA/QAafter/Pion/TPCnclusterPhipi_selected", "TPC ncluster vs phi for selected Pion", kTHnSparseF, {{160, 0, 160, "TPC nCluster"}, {63, 0.0f, 6.28f, "#phi"}});

      //  Mass QA 1D for quick check
      if (cFill1DQAs) {
        histos.add("Result/Data/k892invmass", "Invariant mass of K(892) K^{#pm}#pi^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}});
        histos.add("Result/Data/antik892invmass", "Invariant mass of K(892) K^{#mp}#pi^{#pm}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}});
        histos.add("Result/Data/k892invmassLSPP", "Invariant mass of K(892) Like Sign Method K^{#plus}#pi^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}});   // K+ + Pi
        histos.add("Result/Data/k892invmassLSMM", "Invariant mass of K(892) Like Sign Method K^{#minus}#pi^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}}); // K- + anti-Pi
      }
      // eta phi QA
      if (cFilldeltaEtaPhiPlots) {
        histos.add("QAbefore/deltaEta", "deltaEta of kaon and Pion candidates", HistType::kTH1F, {axisEtaPhi});
        histos.add("QAbefore/deltaPhi", "deltaPhi of kaon and Pion candidates", HistType::kTH1F, {axisEtaPhi});

        histos.add("QAafter/deltaEta", "deltaEta of kaon and Pion candidates", HistType::kTH1F, {axisEtaPhi});
        histos.add("QAafter/deltaPhi", "deltaPhi of kaon and Pion candidates", HistType::kTH1F, {axisEtaPhi});

        histos.add("QAafter/PhiPiafter", "Phi of  Pion candidates", HistType::kTH1F, {axisPhi});
        histos.add("QAafter/PhiKaafter", "Phi of kaon  candidates", HistType::kTH1F, {axisPhi});
      }

      // 3d histogram
      histos.add("Result/Data/h3k892invmass", "Invariant mass of K(892) K^{#pm}#pi^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892});
      histos.add("Result/Data/h3antik892invmass", "Invariant mass of K(892) K^{#mp}#pi^{#pm}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892});
      histos.add("Result/Data/h3k892invmassLSPP", "Invariant mass of K(892) Like Sign Method K^{#plus}#pi^{#plus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892});   // K+ + Pi
      histos.add("Result/Data/h3k892invmassLSMM", "Invariant mass of K(892) Like Sign Method K^{#minus}#pi^{#minus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892}); // K- + anti-Pi
    }

    if (doprocessRotational) {
      if (cFill1DQAs) {
        histos.add("Result/Data/k892InvMassRotation", "Invariant mass of K(892) Like Sign Method K^{#plus}#pi^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}});       // K+ + Pi
        histos.add("Result/Data/antik892InvMassRotation", "Invariant mass of K(892) Like Sign Method K^{#minus}#pi^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}}); // K- + anti-Pi
      }
      histos.add("Result/Data/h3k892InvMassRotation", "Invariant mass of K(892) rotation", kTHnSparseF, {axisMult, axisPt, axisMassk892});
      histos.add("Result/Data/h3antik892InvMassRotation", "Invariant mass of K(892) rotation", kTHnSparseF, {axisMult, axisPt, axisMassk892});
    }
    // Mixed event histograms
    if (doprocessME) {
      if (cFill1DQAs) {
        histos.add("Result/Data/k892invmassME_UnlikeSign", "Invariant mass of K(892) mixed event K^{#pm}#pi^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}});
        histos.add("Result/Data/antik892invmassME_UnlikeSign", "Invariant mass of K(892) mixed event K^{#pm}#pi^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}});
      }
      histos.add("Result/Data/h3k892invmassME_UnlikeSign", "Invariant mass of K(892) mixed event K^{#pm}#pi^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892});
      histos.add("Result/Data/h3antik892invmassME_UnlikeSign", "Invariant mass of K(892) mixed event K^{#pm}#pi^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892});

      if (cFilladditionalMEPlots) {
        if (cFill1DQAs) {
          histos.add("Result/Data/k892invmassME_LSPP", "Invariant mass of K(892) Like Sign Method K^{#plus}#pi^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}});   // K+ + Pi
          histos.add("Result/Data/k892invmassME_LSMM", "Invariant mass of K(892) Like Sign Method K^{#minus}#pi^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassk892}}); // K- + anti-Pi
        }
        histos.add("Result/Data/h3k892invmassME_LSPP", "Invariant mass of K(892) mixed event Like Sign Method K^{#plus}#pi^{#plus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892});   // K+ + Pi
        histos.add("Result/Data/h3k892invmassME_LSMM", "Invariant mass of K(892) mixed event Like Sign Method K^{#minus}#pi^{#minus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassk892}); // K- + anti-Pi
      }
    }

    // MC QA
    if (doprocessEventFactor) {
      histos.add("Event/MultiplicityGenEv", "hMCEventIndices", kTH1D, {axisMult});
      histos.add("Event/MultiplicityRecoEv", "Multiplicity of Reconstructed Events", kTH1D, {axisMult});
    }

    if (doprocessMCGen) {
      histos.add("QA/MC/h2GenEtaPt_beforeanycut", " #eta-#it{p}_{T} distribution of Generated K(892); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_beforeanycut", " #phi-y distribution of Generated K(892); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});
      histos.add("QA/MC/h2GenEtaPt_afterRapcut", " #phi-#it{p}_{T} distribution of Generated K(892); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_afterRapcut", " #phi-y distribution of Generated K(892); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});

      histos.add("Result/MC/Genk892pt", "pT distribution of True MC K(892)0", kTH2F, {axisPt, axisMult});
      histos.add("Result/MC/Genantik892pt", "pT distribution of True MC Anti-K(892)0", kTH2F, {axisPt, axisMult});
    }

    if (doprocessMCRec) {
      histos.add("QA/MC/h2RecoEtaPt_after", " #eta-#it{p}_{T} distribution of Reconstructed K(892); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPt});
      histos.add("QA/MC/h2RecoPhiRapidity_after", " #phi-y distribution of Reconstructed K(892); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});

      histos.add("QA/MC/trkDCAxy_pi", "DCAxy distribution of Pion track candidates", HistType::kTHnSparseF, {axisPt, axisDCAxy});
      histos.add("QA/MC/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTHnSparseF, {axisPt, axisDCAxy});
      histos.add("QA/MC/trkDCAz_pi", "DCAz distribution of Pion track candidates", HistType::kTHnSparseF, {axisPt, axisDCAz});
      histos.add("QA/MC/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTHnSparseF, {axisPt, axisDCAz});
      histos.add("QA/MC/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Pion}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/MC/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Pion}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/MC/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/MC/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});

      histos.add("Result/MC/h3k892Recoinvmass", "Invariant mass of Reconstructed MC K(892)0", kTHnSparseF, {axisMult, axisPt, axisMassk892});
      histos.add("Result/MC/h3antik892Recoinvmass", "Invariant mass of Reconstructed MC Anti-K(892)0", kTHnSparseF, {axisMult, axisPt, axisMassk892});
    }

    if (doprocessSignalLoss) {
      histos.add("Result/SignalLoss/GenTruek892pt_den", "True k892 (den)", kTH2F, {axisPt, axisMult});
      histos.add("Result/SignalLoss/GenTrueantik892pt_den", "True anti-k892 (den)", kTH2F, {axisPt, axisMult});

      histos.add("Result/SignalLoss/GenTruek892pt_num", "True k892 (num)", kTH2F, {axisPt, axisMult});
      histos.add("Result/SignalLoss/GenTrueantik892pt_num", "True anti-k892 (num)", kTH2F, {axisPt, axisMult});
    }

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in Kstar0analysis:";
    histos.print();
  } // end of init

  float massKa = MassKaonCharged;
  float massPi = MassPionCharged;

  // Centralicity estimator selection
  template <typename Coll>
  float centEst(Coll collisions)
  {
    float returnValue = -999.0f;
    switch (centEstimator) {
      case 0:
        returnValue = collisions.centFT0M();
        break;
      case 1:
        returnValue = collisions.centFT0A();
        break;
      case 2:
        returnValue = collisions.centFT0C();
        break;
      default:
        returnValue = collisions.centFT0M();
        break;
    }
    return returnValue;
  }

  template <typename Coll>
  bool isSelected(const Coll& collision, bool fillHist = true)
  {
    auto applyCut = [&](bool enabled, bool condition, int bin) {
      if (!enabled)
        return true;
      if (!condition)
        return false;
      if (fillHist)
        histos.fill(HIST("CollCutCounts"), bin);
      return true;
    };

    if (fillHist)
      histos.fill(HIST("CollCutCounts"), 0);

    if (!applyCut(true, std::abs(collision.posZ()) <= configEvents.cfgEvtZvtx, 1))
      return false;

    if (!applyCut(configEvents.cfgEvtTriggerTVXSel,
                  collision.selection_bit(aod::evsel::kIsTriggerTVX), 2))
      return false;

    if (!applyCut(configEvents.cfgEvtNoTFBorderCut,
                  collision.selection_bit(aod::evsel::kNoTimeFrameBorder), 3))
      return false;

    if (!applyCut(configEvents.cfgEvtNoITSROFrameBorderCut,
                  collision.selection_bit(aod::evsel::kNoITSROFrameBorder), 4))
      return false;

    if (!applyCut(configEvents.cfgEvtIsRCTFlagpassed, rctChecker(collision), 5))
      return false;

    if (!applyCut(configEvents.cfgEvtSel8, collision.sel8(), 6))
      return false;

    if (!applyCut(configEvents.cfgEvtIsINELgt0, collision.isInelGt0(), 7))
      return false;

    if (fillHist)
      histos.fill(HIST("CollCutCounts"), 8);

    return true;
  }

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (configTracks.cDCAr7SigCut && std::abs(track.dcaXY()) > (0.004f + 0.013f / (track.pt()))) // 7 - Sigma cut
      return false;
    if (configTracks.cTPCNClsFound && (track.tpcNClsFound() < configTracks.cMinTPCNClsFound))
      return false;
    if (track.tpcNClsCrossedRows() < configTracks.cfgMinCrossedRows)
      return false;
    if (configTracks.cfgHasTOF && !track.hasTOF())
      return false;
    if (configTracks.cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (configTracks.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (configTracks.cfgPVContributor && !track.isPVContributor())
      return false;
    if (configTracks.cfgGlobalTrack && !track.isGlobalTrack())
      return false;

    return true;
  }

  // LOGF(info, "AFTER: pt: %f, hasTOF: %d, TPCSigma: %f, TOFSigma: %f, boolTPC: %d, boolTOF: %d, bool: %d", pt, candidate.hasTOF(),
  //    candidate.tpcNSigmaPr(), candidate.tofNSigmaPr(), tpcPIDPassed, tofPIDPassed, tpcPIDPassed || tofPIDPassed);

  template <typename T>
  bool ptDependentPidPion(const T& candidate)
  {
    auto vPionTPCPIDpTintv = configPID.pionTPCPIDpTintv.value;
    vPionTPCPIDpTintv.insert(vPionTPCPIDpTintv.begin(), configTracks.cMinPtcut);
    auto vPionTPCPIDcuts = configPID.pionTPCPIDcuts.value;
    auto vPionTOFPIDpTintv = configPID.pionTOFPIDpTintv.value;
    auto vPionTPCTOFCombinedpTintv = configPID.pionTPCTOFCombinedpTintv.value;
    auto vPionTPCTOFCombinedPIDcuts = configPID.pionTPCTOFCombinedPIDcuts.value;
    auto vPionTOFPIDcuts = configPID.pionTOFPIDcuts.value;

    float pt = candidate.pt();
    float ptSwitchToTOF = vPionTPCPIDpTintv.back();
    float tpcNsigmaPi = candidate.tpcNSigmaPi();
    float tofNsigmaPi = candidate.tofNSigmaPi();

    bool tpcPIDPassed = false;

    // TPC PID (interval check)
    for (size_t i = 0; i < vPionTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vPionTPCPIDpTintv[i] && pt < vPionTPCPIDpTintv[i + 1]) {
        if (std::abs(tpcNsigmaPi) < vPionTPCPIDcuts[i])
          tpcPIDPassed = true;
      }
    }

    // TOF bypass option (for QA or MC)
    if (configPID.cByPassTOF) {
      return std::abs(tpcNsigmaPi) < vPionTPCPIDcuts.back();
    }

    // Case 1: No TOF and pt ≤ threshold → accept only via TPC PID
    if (!candidate.hasTOF() && pt <= ptSwitchToTOF) {
      return tpcPIDPassed;
    }

    // Case 2: No TOF but pt > threshold → reject
    if (!candidate.hasTOF() && pt > ptSwitchToTOF) {
      return false;
    }

    // Case 3: Has TOF → use TPC + TOF (square or circular)
    if (candidate.hasTOF()) {
      if (configPID.cPIDcutType == SquareType) {
        // Rectangular cut
        for (size_t i = 0; i < vPionTOFPIDpTintv.size(); ++i) {
          if (pt < vPionTOFPIDpTintv[i]) {
            if (std::abs(tofNsigmaPi) < vPionTOFPIDcuts[i] &&
                std::abs(tpcNsigmaPi) < vPionTPCPIDcuts.back())
              return true;
          }
        }
      } else if (configPID.cPIDcutType == CircularType) {
        // Circular cut
        for (size_t i = 0; i < vPionTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vPionTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 =
              tpcNsigmaPi * tpcNsigmaPi +
              tofNsigmaPi * tofNsigmaPi;
            if (combinedSigma2 < vPionTPCTOFCombinedPIDcuts[i] * vPionTPCTOFCombinedPIDcuts[i])
              return true;
          }
        }
      }
    }

    return false;
  }

  template <typename T>
  bool ptDependentPidKaon(const T& candidate)
  {
    auto vKaonTPCPIDpTintv = configPID.kaonTPCPIDpTintv.value;
    vKaonTPCPIDpTintv.insert(vKaonTPCPIDpTintv.begin(), configTracks.cMinPtcut);
    auto vKaonTPCPIDcuts = configPID.kaonTPCPIDcuts.value;
    auto vKaonTOFPIDpTintv = configPID.kaonTOFPIDpTintv.value;
    auto vKaonTPCTOFCombinedpTintv = configPID.kaonTPCTOFCombinedpTintv.value;
    auto vKaonTPCTOFCombinedPIDcuts = configPID.kaonTPCTOFCombinedPIDcuts.value;
    auto vKaonTOFPIDcuts = configPID.kaonTOFPIDcuts.value;

    float pt = candidate.pt();
    float ptSwitchToTOF = vKaonTPCPIDpTintv.back();
    float tpcNsigmaKa = candidate.tpcNSigmaKa();
    float tofNsigmaKa = candidate.tofNSigmaKa();

    bool tpcPIDPassed = false;

    // TPC PID interval-based check
    for (size_t i = 0; i < vKaonTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vKaonTPCPIDpTintv[i] && pt < vKaonTPCPIDpTintv[i + 1]) {
        if (std::abs(tpcNsigmaKa) < vKaonTPCPIDcuts[i]) {
          tpcPIDPassed = true;
          break;
        }
      }
    }

    // TOF bypass option
    if (configPID.cByPassTOF) {
      return std::abs(tpcNsigmaKa) < vKaonTPCPIDcuts.back();
    }

    // Case 1: No TOF and pt ≤ ptSwitch → use TPC-only
    if (!candidate.hasTOF() && pt <= ptSwitchToTOF) {
      return tpcPIDPassed;
    }

    // Case 2: No TOF but pt > ptSwitch → reject
    if (!candidate.hasTOF() && pt > ptSwitchToTOF) {
      return false;
    }

    // Case 3: TOF is available → apply TPC+TOF PID logic
    if (candidate.hasTOF()) {
      if (configPID.cPIDcutType == SquareType) {
        // Rectangular cut
        for (size_t i = 0; i < vKaonTOFPIDpTintv.size(); ++i) {
          if (pt < vKaonTOFPIDpTintv[i]) {
            if (std::abs(tofNsigmaKa) < vKaonTOFPIDcuts[i] &&
                std::abs(tpcNsigmaKa) < vKaonTPCPIDcuts.back()) {
              return true;
            }
          }
        }
      } else if (configPID.cPIDcutType == CircularType) {
        // Circular cut
        for (size_t i = 0; i < vKaonTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vKaonTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 = tpcNsigmaKa * tpcNsigmaKa +
                                   tofNsigmaKa * tofNsigmaKa;
            if (combinedSigma2 < vKaonTPCTOFCombinedPIDcuts[i] * vKaonTPCTOFCombinedPIDcuts[i]) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  auto static constexpr MinPtforProtonRejection = 1.0f;
  auto static constexpr MaxPtforProtonRejection = 2.0f;
  auto static constexpr MaxnSigmaforProtonRejection = 2.0f;

  template <typename T>
  bool rejectProton(const T& candidate)
  {
    if (candidate.pt() > MinPtforProtonRejection && candidate.pt() < MaxPtforProtonRejection && !candidate.hasTOF() && candidate.tpcNSigmaPr() < MaxnSigmaforProtonRejection) {
      return false;
    }
    return true;
  }

  auto static constexpr MaxNok892Daughters = 2;

  template <bool IsData, bool IsRot, bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto centrality = centEst(collision);

    // Multiplicity correlation calibration plots
    if (cFillMultQA) {
      if constexpr (IsData) {
        histos.fill(HIST("MultCalib/centGloPVpr"), centrality, dTracks1.size(), collision.multNTracksPV());
        histos.fill(HIST("MultCalib/centGloPVka"), centrality, dTracks2.size(), collision.multNTracksPV());
      }
    }

    if (cFilladditionalQAeventPlots) {
      if constexpr (IsData) {
        histos.fill(HIST("QAevent/hnTrksSameE"), dTracks1.size());
      } else if constexpr (IsMix) {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), centrality);
        histos.fill(HIST("QAevent/hnTrksMixedE"), dTracks1.size());
      }
    }
    // LOG(info) << "After pass, Collision index:" << collision.index() << "multiplicity: " << collision.centFT0M() << endl;

    LorentzVectorPtEtaPhiMass lDecayDaughter1, lDecayDaughter2, lResonance, ldaughterRot, lResonanceRot;

    for (const auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      //// Initialize variables
      // Trk1: Pion
      auto isTrk1hasTOF = trk1.hasTOF();
      auto trk1ptPi = trk1.pt();
      auto trk1etaPi = trk1.eta();
      auto trk1phiPi = trk1.phi();
      auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.0f;

      // Trk2: Kaon
      auto isTrk2hasTOF = trk2.hasTOF();
      auto trk2ptKa = trk2.pt();
      auto trk2etaKa = trk2.eta();
      auto trk2phiKa = trk2.phi();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.0f;

      auto deltaEta = 0.0f;
      auto deltaPhi = 0.0f;

      if (cfgUseDeltaEtaPhiCuts) {
        deltaEta = std::abs(trk1etaPi - trk2etaKa);
        deltaPhi = std::abs(trk1phiPi - trk2phiKa);
        deltaPhi = (deltaPhi > o2::constants::math::PI) ? (o2::constants::math::TwoPI - deltaPhi) : deltaPhi;
        if (deltaEta >= cMaxDeltaEtaCut)
          continue;
        if (deltaPhi >= cMaxDeltaPhiCut)
          continue;
      }

      //// QA plots before the selection
      //  --- Track QA all
      if constexpr (IsData) {
        if (cFillTrackQA) {
          histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_pi_all"), centrality, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_pi_all"), centrality, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_pi_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
          if (!isTrk1hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_pi_all"), trk1ptPi, trk1NSigmaPiTPC);
          }
          histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTOF);
            histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
          }
          if (!isTrk2hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_ka_all"), trk2ptKa, trk2NSigmaKaTPC);
          }

          histos.fill(HIST("QA/QAbefore/Track/dcaZ_all"), trk1ptPi, trk1.dcaZ());
          histos.fill(HIST("QA/QAbefore/Track/dcaXY_all"), trk1ptPi, trk1.dcaXY());
          histos.fill(HIST("QA/QAbefore/Track/TPC_CR_all"), trk1ptPi, trk1.tpcNClsCrossedRows());
          histos.fill(HIST("QA/QAbefore/Track/pT_all"), trk1ptPi);
          histos.fill(HIST("QA/QAbefore/Track/eta_all"), trk1etaPi);
        }
        if (cFilldeltaEtaPhiPlots) {
          histos.fill(HIST("QAbefore/deltaEta"), deltaEta);
          histos.fill(HIST("QAbefore/deltaPhi"), deltaPhi);
        }
      }

      //// Apply the pid selection
      if (crejectProton && rejectProton(trk2)) // to remove proton contamination from the kaon track
        continue;

      if (!ptDependentPidPion(trk1) || !ptDependentPidKaon(trk2))
        continue;

      //// QA plots after the selection
      if constexpr (IsData) { //  --- PID QA Pion
        histos.fill(HIST("QA/QAafter/Pion/TPC_Nsigma_pi_selected"), centrality, trk1ptPi, trk1NSigmaPiTPC);
        histos.fill(HIST("QA/QAafter/Pion/TPC_Signal_pi_selected"), trk1.tpcInnerParam(), trk1.tpcSignal());
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Pion/TOF_Nsigma_pi_selected"), centrality, trk1ptPi, trk1NSigmaPiTOF);
          histos.fill(HIST("QA/QAafter/Pion/TOF_TPC_Map_pi_selected"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Pion/TPC_Nsigma_pi_TPConly_selected"), trk1ptPi, trk1NSigmaPiTPC);
        }
        histos.fill(HIST("QA/QAafter/Pion/dcaZ_selected"), trk1ptPi, trk1.dcaZ());
        histos.fill(HIST("QA/QAafter/Pion/dcaXY_selected"), trk1ptPi, trk1.dcaXY());
        histos.fill(HIST("QA/QAafter/Pion/TPC_CR_selected"), trk1ptPi, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Pion/pT_selected"), trk1ptPi);
        histos.fill(HIST("QA/QAafter/Pion/eta_selected"), trk1etaPi);
        histos.fill(HIST("QA/QAafter/Pion/TPCnclusterPhipi_selected"), trk1.tpcNClsFound(), trk1phiPi);

        //  --- PID QA Kaon
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_selected"), centrality, trk2ptKa, trk2NSigmaKaTPC);
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Signal_ka_selected"), trk2.tpcInnerParam(), trk2.tpcSignal());
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TOF_Nsigma_ka_selected"), centrality, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAafter/Kaon/TOF_TPC_Map_ka_selected"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly_selected"), trk2ptKa, trk2NSigmaKaTPC);
        }
        histos.fill(HIST("QA/QAafter/Kaon/dcaZ_selected"), trk2ptKa, trk2.dcaZ());
        histos.fill(HIST("QA/QAafter/Kaon/dcaXY_selected"), trk2ptKa, trk2.dcaXY());
        histos.fill(HIST("QA/QAafter/Kaon/TPC_CR_selected"), trk2ptKa, trk2.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Kaon/pT_selected"), trk2ptKa);
        histos.fill(HIST("QA/QAafter/Kaon/eta_selected"), trk2etaKa);
        histos.fill(HIST("QA/QAafter/Kaon/TPCnclusterPhika_selected"), trk2.tpcNClsFound(), trk2phiKa);

        if (cFilldeltaEtaPhiPlots) {
          histos.fill(HIST("QAafter/PhiPrafter"), trk1phiPi);
          histos.fill(HIST("QAafter/PhiKaafter"), trk2phiKa);
          histos.fill(HIST("QAafter/deltaEta"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhi"), deltaPhi);
        }
      }

      //// Resonance reconstruction
      lDecayDaughter1 = LorentzVectorPtEtaPhiMass(trk1ptPi, trk1etaPi, trk1phiPi, massPi);
      lDecayDaughter2 = LorentzVectorPtEtaPhiMass(trk2ptKa, trk2etaKa, trk2phiKa, massKa);

      // Apply kinematic opening angle cut
      if (cUseOpeningAngleCut) {
        auto v1 = lDecayDaughter1.Vect();
        auto v2 = lDecayDaughter2.Vect();
        float alpha = std::acos(v1.Dot(v2) / (v1.R() * v2.R()));
        if (alpha > cMinOpeningAngle && alpha < cMaxOpeningAngle)
          continue;
      }

      lResonance = lDecayDaughter1 + lDecayDaughter2;

      auto resonanceMass = lResonance.M();
      auto resonancePt = lResonance.Pt();
      auto resonanceRapidity = lResonance.Rapidity();

      if constexpr (IsData || IsMix) {
        // Rapidity cut
        if (std::abs(resonanceRapidity) > configTracks.cfgCutRapidity)
          continue;
      }

      if (cfgUseCutsOnMother) {
        if (resonancePt >= cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (resonanceMass >= cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      }

      if (cFilladditionalQAeventPlots) {
        if constexpr (IsData) {
          histos.fill(HIST("QAevent/hPairsCounterSameE"), 1.0f);
        } else if (IsMix) {
          histos.fill(HIST("QAevent/hPairsCounterMixedE"), 1.0f);
        }
      }

      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (IsRot) {
          for (int i = 0; i < configBkg.cNofRotations; i++) {
            float theta = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / configBkg.rotationalcut, o2::constants::math::PI + o2::constants::math::PI / configBkg.rotationalcut);
            if (configBkg.cfgRotPi) {
              ldaughterRot = LorentzVectorPtEtaPhiMass(trk1ptPi, trk1etaPi, trk1phiPi + theta, massPi);
              lResonanceRot = ldaughterRot + lDecayDaughter2;
            } else {
              ldaughterRot = LorentzVectorPtEtaPhiMass(trk2ptKa, trk2etaKa, trk2phiKa + theta, massKa);
              lResonanceRot = lDecayDaughter1 + ldaughterRot;
            }
            auto resonanceRotMass = lResonanceRot.M();
            auto resonanceRotPt = lResonanceRot.Pt();

            // Rapidity cut
            if (std::abs(lResonanceRot.Rapidity()) >= configTracks.cfgCutRapidity)
              continue;

            if (cfgUseCutsOnMother) {
              if (resonanceRotPt >= cMaxPtMotherCut) // excluding candidates in overflow
                continue;
              if (resonanceRotMass >= cMaxMinvMotherCut) // excluding candidates in overflow
                continue;
            }
            if (trk1.sign() < 0) {
              if (cFill1DQAs) {
                histos.fill(HIST("Result/Data/k892InvMassRotation"), resonanceRotMass);
              }
              histos.fill(HIST("Result/Data/h3k892InvMassRotation"), centrality, resonanceRotPt, resonanceRotMass);
            } else if (trk1.sign() > 0) {
              if (cFill1DQAs) {
                histos.fill(HIST("Result/Data/antik892InvMassRotation"), resonanceRotMass);
              }
              histos.fill(HIST("Result/Data/h3antik892InvMassRotation"), centrality, resonanceRotPt, resonanceRotMass);
            }
          }
        }
        if constexpr (IsData) {
          if (trk1.sign() < 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/k892invmass"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3k892invmass"), centrality, resonancePt, resonanceMass);
          } else if (trk1.sign() > 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/antik892invmass"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3antik892invmass"), centrality, resonancePt, resonanceMass);
          }
        } else if (IsMix) {
          if (trk1.sign() < 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/k892invmassME_UnlikeSign"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3k892invmassME_UnlikeSign"), centrality, resonancePt, resonanceMass);
          } else if (trk1.sign() > 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/antik892invmassME_UnlikeSign"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3antik892invmassME_UnlikeSign"), centrality, resonancePt, resonanceMass);
          }
        }

        // MC
        if constexpr (IsMC) {
          // now we do mc true
          // ------ Temporal lambda function to prevent error in build
          auto getMothersIndeces = [&](auto const& theMcParticle) {
            std::vector<int> lMothersIndeces{};
            for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
              LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
              lMothersIndeces.push_back(lMother.globalIndex());
            }
            return lMothersIndeces;
          };
          auto getMothersPDGCodes = [&](auto const& theMcParticle) {
            std::vector<int> lMothersPDGs{};
            for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
              LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
              lMothersPDGs.push_back(lMother.pdgCode());
            }
            return lMothersPDGs;
          };
          // ------
          std::vector<int> motherstrk1 = {-1, -1};
          std::vector<int> mothersPDGtrk1 = {-1, -1};

          std::vector<int> motherstrk2 = {-1, -1};
          std::vector<int> mothersPDGtrk2 = {-1, -1};

          //
          // Get the MC particle
          const auto& mctrk1 = trk1.mcParticle();
          if (mctrk1.has_mothers()) {
            motherstrk1 = getMothersIndeces(mctrk1);
            mothersPDGtrk1 = getMothersPDGCodes(mctrk1);
          }
          while (motherstrk1.size() > MaxNok892Daughters) {
            motherstrk1.pop_back();
            mothersPDGtrk1.pop_back();
          }

          const auto& mctrk2 = trk2.mcParticle();
          if (mctrk2.has_mothers()) {
            motherstrk2 = getMothersIndeces(mctrk2);
            mothersPDGtrk2 = getMothersPDGCodes(mctrk2);
          }
          while (motherstrk2.size() > MaxNok892Daughters) {
            motherstrk2.pop_back();
            mothersPDGtrk2.pop_back();
          }

          if (std::abs(mctrk1.pdgCode()) != kPiPlus || std::abs(mctrk2.pdgCode()) != kKPlus)
            continue;

          if (motherstrk1[0] != motherstrk2[0]) // Same mother
            continue;

          if (std::abs(mothersPDGtrk1[0]) != Pdg::kK0Star892)
            continue;

          // LOGF(info, "mother trk1 id: %d, mother trk1: %d, trk1 id: %d, trk1 pdgcode: %d, mother trk2 id: %d, mother trk2: %d, trk2 id: %d, trk2 pdgcode: %d", motherstrk1[0], mothersPDGtrk1[0], trk1.globalIndex(), mctrk1.pdgCode(), motherstrk2[0], mothersPDGtrk2[0], trk2.globalIndex(), mctrk2.pdgCode());

          if (cUseRapcutMC && std::abs(resonanceRapidity) > configTracks.cfgCutRapidity) // rapidity cut
            continue;

          histos.fill(HIST("QA/MC/h2RecoEtaPt_after"), lResonance.Eta(), resonancePt);
          histos.fill(HIST("QA/MC/h2RecoPhiRapidity_after"), lResonance.Phi(), resonanceRapidity);

          // Track selection check.
          histos.fill(HIST("QA/MC/trkDCAxy_pi"), trk1ptPi, trk1.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAxy_ka"), trk2ptKa, trk2.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAz_pi"), trk1ptPi, trk1.dcaZ());
          histos.fill(HIST("QA/MC/trkDCAz_ka"), trk2ptKa, trk2.dcaZ());

          histos.fill(HIST("QA/MC/TPC_Nsigma_pi_all"), centrality, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_pi_all"), centrality, trk1ptPi, trk1NSigmaPiTOF);
          }
          histos.fill(HIST("QA/MC/TPC_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTOF);
          }

          // MC histograms
          if (mothersPDGtrk1[0] > 0) {
            histos.fill(HIST("Result/MC/h3k892Recoinvmass"), centrality, resonancePt, resonanceMass);
          } else {
            histos.fill(HIST("Result/MC/h3antik892Recoinvmass"), centrality, resonancePt, resonanceMass);
          }
        }
      } else {
        if constexpr (IsData) {
          // Like sign pair ++
          if (trk1.sign() > 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/k892invmassLSPP"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3k892invmassLSPP"), centrality, resonancePt, resonanceMass);
          } else { // Like sign pair --
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/k892invmassLSMM"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3k892invmassLSMM"), centrality, resonancePt, resonanceMass);
          }
        } else if (IsMix) {
          if (cFilladditionalMEPlots) {
            // Like sign pair ++
            if (trk1.sign() > 0) {
              if (cFill1DQAs) {
                histos.fill(HIST("Result/Data/k892invmassME_LSPP"), resonanceMass);
              }
              histos.fill(HIST("Result/Data/h3k892invmassME_LSPP"), centrality, resonancePt, resonanceMass);
            } else { // Like sign pair --
              if (cFill1DQAs) {
                histos.fill(HIST("Result/Data/k892invmassME_LSMM"), resonanceMass);
              }
              histos.fill(HIST("Result/Data/h3k892invmassME_LSMM"), centrality, resonancePt, resonanceMass);
            }
          }
        }
      }
    }
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks)
  {
    if (!isSelected(collision)) // Default event selection
      return;

    auto centrality = centEst(collision);

    histos.fill(HIST("Event/posZ"), collision.posZ());
    histos.fill(HIST("Event/centFT0M"), centrality);

    fillHistograms<true, false, false, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Kstar0analysis, processData, "Process Event for data without partition", false);

  void processRotational(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    if (!isSelected(collision, false)) // Default event selection
      return;

    if (!collision.isInelGt0()) // <--
      return;

    fillHistograms<false, true, false, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Kstar0analysis, processRotational, "Process Rotational Background", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<collision::PosZ, cent::CentFT0M>;

  void processME(EventCandidates const& collision,
                 TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);

    BinningTypeVtxZT0M colBinning{{configBkg.cfgVtxBins, configBkg.cfgMultPercentileBins}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVtxZT0M> pairs{colBinning, configBkg.nEvtMixing, -1, collision, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      // LOGF(info, "Mixed event collisions: (%d, %d)", collision1.globalIndex(), collision2.globalIndex());

      // for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
      //  LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.index(), t2.index(), collision1.index(), collision2.index());
      //  }

      if (!isSelected(collision1, false)) // Default event selection
        continue;

      if (!isSelected(collision2, false)) // Default event selection
        continue;

      if (!collision1.isInelGt0()) // <--
        continue;

      if (!collision2.isInelGt0()) // <--
        continue;

      if (cFilladditionalQAeventPlots) {
        // Fill histograms for the characteristics of the *mixed* events (collision1 and collision2)
        // This will show the distribution of events that are actually being mixed.
        if (cFill1DQAs) {
          histos.fill(HIST("QAevent/hMixPool_VtxZ"), collision1.posZ());
          histos.fill(HIST("QAevent/hMixPool_Multiplicity"), collision1.centFT0M());
        }
        histos.fill(HIST("QAevent/hMixPool_VtxZ_vs_Multiplicity"), collision1.posZ(), collision1.centFT0M());

        // You might also want to fill for collision2 if you want to see both partners' distributions
        // histos.fill(HIST("QAevent/hMixPool_VtxZ"), collision2.posZ());
        // histos.fill(HIST("QAevent/hMixPool_Multiplicity"), collision2.centFT0M());
        // histos.fill(HIST("QAevent/hMixPool_VtxZ_vs_Multiplicity"), collision2.posZ(), collision2.centFT0M());
      }
      fillHistograms<false, false, false, true>(collision1, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(Kstar0analysis, processME, "Process EventMixing light without partition", false);

  void processMCRec(MCEventCandidates::iterator const& collision,
                    aod::McCollisions const&,
                    MCTrackCandidates const& tracks, aod::McParticles const&)
  {
    if (!collision.has_mcCollision())
      return;

    if (!isSelected(collision))
      return;

    auto centrality = centEst(collision);

    histos.fill(HIST("Event/posZ"), collision.posZ());
    histos.fill(HIST("Event/centFT0M"), centrality);

    fillHistograms<false, false, true, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Kstar0analysis, processMCRec, "Process Event for MC Rec without partition", false);

  Partition<aod::McParticles> selectedMCParticles = (nabs(aod::mcparticle::pdgCode) == static_cast<int>(Pdg::kK0Star892)); // K(892)

  void processMCGen(MCEventCandidates::iterator const& collision, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    if (!collision.has_mcCollision())
      return;

    bool isInAfterAllCuts = isSelected(collision, false);

    auto mcPartsAll = mcParticles.sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    // bool isTrueINELgt0 = pwglf::isINELgt0mc(mcPartsAll, pdg);
    bool isTrueINELgt0 = collision.isInelGt0(); // <--

    auto centrality = centEst(collision);

    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());

    // Not related to the real collisions
    for (const auto& part : mcParts) { // loop over all K(892) particles

      std::vector<int> daughterPDGs;
      if (part.has_daughters()) {
        auto daughter01 = mcParticles.rawIteratorAt(part.daughtersIds()[0] - mcParticles.offset());
        auto daughter02 = mcParticles.rawIteratorAt(part.daughtersIds()[1] - mcParticles.offset());
        daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
      } else {
        daughterPDGs = {-1, -1};
      }

      bool pass1 = std::abs(daughterPDGs[0]) == kKPlus || std::abs(daughterPDGs[1]) == kKPlus;   // At least one decay to Kaon
      bool pass2 = std::abs(daughterPDGs[0]) == kPiPlus || std::abs(daughterPDGs[1]) == kPiPlus; // At least one decay to Pion

      // Checking if we have both decay products
      if (!pass1 || !pass2)
        continue;

      // LOGF(info, "Part PDG: %d", part.pdgCode(), "DAU_ID1: %d", pass1, "DAU_ID2: %d", pass2);

      histos.fill(HIST("QA/MC/h2GenEtaPt_beforeanycut"), part.eta(), part.pt());
      histos.fill(HIST("QA/MC/h2GenPhiRapidity_beforeanycut"), part.phi(), part.y());

      if (cUseRapcutMC && std::abs(part.y()) > configTracks.cfgCutRapidity) // rapidity cut
        continue;

      if (cfgUseDaughterEtaCutMC) {
        for (auto const& daughters : part.daughters_as<aod::McParticles>()) {
          if (std::fabs(daughters.eta()) > configTracks.cfgCutEta)
            continue; // eta cut for daughters
        } // loop over daughters
      }

      histos.fill(HIST("QA/MC/h2GenEtaPt_afterRapcut"), part.eta(), part.pt());
      histos.fill(HIST("QA/MC/h2GenPhiRapidity_afterRapcut"), part.phi(), part.y());

      if (isInAfterAllCuts && isTrueINELgt0) // after all event selection && INEL>0
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genk892pt"), part.pt(), centrality);
        else
          histos.fill(HIST("Result/MC/Genantik892pt"), part.pt(), centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstar0analysis, processMCGen, "Process Event for MC only", false);

  void processEventFactor(MCEventCandidates const& collisions, soa::Join<aod::McCollisions, aod::McCentFT0Ms> const& mcCollisions, aod::McParticles const& mcParticles)
  {
    // Loop on reconstructed collisions
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cacheMC);

      bool isTrueINELgt0 = pwglf::isINELgt0mc(particlesInCollision, pdg);
      bool isInAfterAllCuts = isSelected(collision, false);

      float centrality = mcCollision.centFT0M();

      if (isTrueINELgt0 && isInAfterAllCuts)
        histos.fill(HIST("Event/MultiplicityRecoEv"), centrality);
    }

    // Loop on generated collisions to fill the event factor for the INEL>0 correction
    for (const auto& mccolls : mcCollisions) {
      float centrality = mccolls.centFT0M();
      bool inVtx10 = std::abs(mccolls.posZ()) <= configEvents.cfgEvtZvtx;

      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolls.globalIndex(), cacheMC);
      bool isTrueINELgt0 = pwglf::isINELgt0mc(particlesInCollision, pdg); // QA for Trigger efficiency

      if (inVtx10 && isTrueINELgt0)
        histos.fill(HIST("Event/MultiplicityGenEv"), centrality);
    }
  }
  PROCESS_SWITCH(Kstar0analysis, processEventFactor, "Process Event factor", false);

  void processSignalLoss(MCEventCandidates const& collisions, soa::Join<aod::McCollisions, aod::McCentFT0Ms> const& mcCollisions, aod::McParticles const& mcParticles)
  {
    // Loop on reconstructed collisions
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cacheMC);

      bool isTrueINELgt0 = pwglf::isINELgt0mc(particlesInCollision, pdg);
      bool isInAfterAllCuts = isSelected(collision, false);
      bool inVtx10 = std::abs(mcCollision.posZ()) <= configEvents.cfgEvtZvtx;

      float centrality = mcCollision.centFT0M();

      // ===== NUM =====
      if (!(inVtx10 && isTrueINELgt0))
        continue;

      if (!isInAfterAllCuts)
        continue;

      for (const auto& part : particlesInCollision) {

        if (cUseRapcutMC && std::abs(part.y()) > configTracks.cfgCutRapidity)
          continue;

        float pt = part.pt();

        // kaon
        if (std::abs(part.pdgCode()) == Pdg::kK0Star892) {
          if (part.pdgCode() > 0)
            histos.fill(HIST("Result/SignalLoss/GenTruek892pt_num"), pt, centrality);
          else
            histos.fill(HIST("Result/SignalLoss/GenTrueantik892pt_num"), pt, centrality);
        }
      }
    }

    // Loop on generated collisions to fill the event factor for the INEL>0 correction
    for (const auto& mccolls : mcCollisions) {
      float centrality = mccolls.centFT0M();

      bool inVtx10 = std::abs(mccolls.posZ()) <= configEvents.cfgEvtZvtx;

      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolls.globalIndex(), cacheMC);
      bool isTrueINELgt0 = pwglf::isINELgt0mc(particlesInCollision, pdg);

      if (!(inVtx10 && isTrueINELgt0))
        continue;

      for (const auto& part : particlesInCollision) {

        if (cUseRapcutMC && std::abs(part.y()) > configTracks.cfgCutRapidity)
          continue;

        // =========================
        // ===== K(892) ======
        // =========================
        if (std::abs(part.pdgCode()) == Pdg::kK0Star892) {

          std::vector<int> daughterPDGs;
          if (part.has_daughters()) {
            auto daughter01 = mcParticles.rawIteratorAt(part.daughtersIds()[0] - mcParticles.offset());
            auto daughter02 = mcParticles.rawIteratorAt(part.daughtersIds()[1] - mcParticles.offset());
            daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
          } else {
            daughterPDGs = {-1, -1};
          }

          bool pass1 = std::abs(daughterPDGs[0]) == kKPlus || std::abs(daughterPDGs[1]) == kKPlus;   // At least one decay to Kaon
          bool pass2 = std::abs(daughterPDGs[0]) == kPiPlus || std::abs(daughterPDGs[1]) == kPiPlus; // At least one decay to Pion

          // Checking if we have both decay products
          if (!pass1 || !pass2)
            continue;

          if (part.pdgCode() > 0)
            histos.fill(HIST("Result/SignalLoss/GenTruek892pt_den"), part.pt(), centrality);
          else
            histos.fill(HIST("Result/SignalLoss/GenTrueantik892pt_den"), part.pt(), centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstar0analysis, processSignalLoss, "Process SignalLoss", false);
}; // Kstar0analysis
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Kstar0analysis>(cfgc)};
};
