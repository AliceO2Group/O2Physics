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

/// \file lambda1520analysisinpp.cxx
/// \brief This standalone task reconstructs track-track decay of lambda(1520) resonance candidate
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

#include "PWGLF/Utils/collisionCuts.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TPDGCode.h"
#include "TRandom.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using LorentzVectorPtEtaPhiMass = ROOT::Math::PtEtaPhiMVector;

enum {
  Inel = 1,
  Inel10,
  Inelg0,
  Inelg010,
  Trig,
  Trig10,
  TrigINELg0,
  TrigINELg010,
  Sel8,
  Sel810,
  Sel8INELg0,
  Sel8INELg010,
  AllCuts,
  AllCuts10,
  AllCutsINELg0,
  AllCutsINELg010,
};

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

struct Lambda1520analysisinpp {
  // Define slice per Resocollision
  SliceCache cache;
  Preslice<Tracks> perCollision = o2::aod::track::collisionId;
  Preslice<McParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<framework::O2DatabasePDG> pdg;

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;

  struct : ConfigurableGroup {
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.0f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<int> cfgEvtOccupancyInTimeRangeMax{"cfgEvtOccupancyInTimeRangeMax", -1, "Evt sel: maximum track occupancy"};
    Configurable<int> cfgEvtOccupancyInTimeRangeMin{"cfgEvtOccupancyInTimeRangeMin", -1, "Evt sel: minimum track occupancy"};
    Configurable<bool> cfgEvtSel8{"cfgEvtSel8", false, "Evt Sel 8 check for offline selection"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", true, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtNoTFBorderCut{"cfgEvtNoTFBorderCut", true, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtIsVertexITSTPC{"cfgEvtIsVertexITSTPC", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
    Configurable<bool> cfgEvtIsGoodZvtxFT0vsPV{"cfgEvtIsGoodZvtxFT0vsPV", true, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> cfgEvtNoSameBunchPileup{"cfgEvtNoSameBunchPileup", false, "Evt sel: apply pileup rejection"};
    Configurable<bool> cfgEvtNoITSROFrameBorderCut{"cfgEvtNoITSROFrameBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtNoCollInTimeRangeStandard{"cfgEvtNoCollInTimeRangeStandard", false, "Evt sel: apply NoNoCollInTimeRangeStandard"};
    Configurable<bool> cfgEvtIsVertexTOFmatched{"cfgEvtIsVertexTOFmatched", true, "kIsVertexTOFmatched: apply vertex TOF matched"};
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

    // Proton
    Configurable<std::vector<float>> protonTPCPIDpTintv{"protonTPCPIDpTintv", {0.9f}, "pT intervals for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> protonTPCPIDcuts{"protonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> protonTOFPIDpTintv{"protonTOFPIDpTintv", {999.0f}, "pT intervals for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> protonTOFPIDcuts{"protonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> protonTPCTOFCombinedpTintv{"protonTPCTOFCombinedpTintv", {999.0f}, "pT intervals for Proton TPC-TOF PID cuts"};
    Configurable<std::vector<float>> protonTPCTOFCombinedPIDcuts{"protonTPCTOFCombinedPIDcuts", {2}, "nSigma list for Proton TPC-TOF PID cuts"};
  } configPID;

  struct : ConfigurableGroup {
    /// Event Mixing
    Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgMultPercentileBins{"cfgMultPercentileBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - multiplicity"};

    // Rotational background
    Configurable<int> rotationalcut{"rotationalcut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
    Configurable<int> cNofRotations{"cNofRotations", 3, "Number of random rotations in the rotational background"};
    Configurable<bool> cfgRotPr{"cfgRotPr", true, "rotate Proton"};
  } configBkg;

  // Additional purity check
  Configurable<bool> crejectPion{"crejectPion", false, "Switch to turn on/off pion contamination"};
  Configurable<bool> cUseOpeningAngleCut{"cUseOpeningAngleCut", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<float> cMinOpeningAngle{"cMinOpeningAngle", 1.4f, "Minimum opening angle between daughters"};
  Configurable<float> cMaxOpeningAngle{"cMaxOpeningAngle", 2.4f, "Maximum opening angle between daughters"};
  Configurable<bool> cfgUseDeltaEtaPhiCuts{"cfgUseDeltaEtaPhiCuts", false, "Switch to turn on/off delta eta and delta phi cuts"};
  Configurable<bool> cfgUseDaughterEtaCutMC{"cfgUseDaughterEtaCutMC", false, "Switch to turn on/off eta cuts for daughters in MC"};

  // MC selection cut
  Configurable<float> cEtacutMC{"cEtacutMC", 0.5f, "MC eta cut"};
  Configurable<bool> cUseRapcutMC{"cUseRapcutMC", true, "MC eta cut"};
  Configurable<bool> cUseEtacutMC{"cUseEtacutMC", true, "MC eta cut"};

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
  // Filter centralityFilter = nabs(aod::cent::centFT0C) <= cfg_Event_CentralityMax;
  // Filter triggerFilter = (o2::aod::evsel::sel8 == true);

  Filter acceptanceFilter = (nabs(aod::track::eta) < configTracks.cfgCutEta && nabs(aod::track::pt) > configTracks.cMinPtcut) &&
                            (nabs(aod::track::dcaXY) < configTracks.cMaxDCArToPVcut) && (nabs(aod::track::dcaZ) < configTracks.cMaxDCAzToPVcut);

  // Filter tofPIDFilter = aod::track::tofExpMom < 0.0f || ((aod::track::tofExpMom > 0.0f) && ( (nabs(aod::pidtof::tofNSigmaPi) < configPID.pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaKa) < configPID.pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaPr) < configPID.pidnSigmaPreSelectionCut))); // TOF
  // Filter tpcPIDFilter = nabs(aod::pidtpc::tpcNSigmaPi) < configPID.pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaKa) < configPID.pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaPr) < configPID.pidnSigmaPreSelectionCut;                                                                             // TPC
  /*  Filter trackFilter = (configTracks.trackSelection == AllTracks) ||
                        ((configTracks.trackSelection == GlobalTracks) && requireGlobalTrackInFilter()) ||
                        ((configTracks.trackSelection == GlobalTracksWoPtEta) && requireGlobalTrackWoPtEtaInFilter()) ||
                        ((configTracks.trackSelection == GlobalTracksWoDCA) && requireGlobalTrackWoDCAInFilter()) ||
                        ((configTracks.trackSelection == QualityTracks) && requireQualityTracksInFilter()) ||
                        ((configTracks.trackSelection == InAcceptanceTracks) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));
  */
  // Filter primarytrackFilter = requirePVContributor() && requirePrimaryTrack() && requireGlobalTrackWoDCA();

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>>;
  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Filtered<soa::Join<TrackCandidates, aod::McTrackLabels>>;

  /// Figures
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.12f, 0.14f, 0.16f, 0.18f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.45f, 0.5f, 0.55f, 0.6f, 0.65f, 0.7f, 0.75f, 0.8f, 0.85f, 0.9f, 0.95f, 1.0f, 1.1f, 1.2f, 1.25f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.75f, 1.8f, 1.9f, 2.0f, 2.1f, 2.2f, 2.3f, 2.4f, 2.5f, 2.6f, 2.7f, 2.8f, 2.9f, 3.0f, 3.1f, 3.2f, 3.3f, 3.4f, 3.6f, 3.7f, 3.8f, 3.9f, 4.0f, 4.1f, 4.2f, 4.5f, 4.6f, 4.8f, 4.9f, 5.0f, 5.5f, 5.6f, 6.0f, 6.4f, 6.5f, 7.0f, 7.2f, 8.0f, 9.0f, 9.5f, 9.6f, 10.0f, 11.0f, 11.5f, 12.0f, 13.0f, 14.0f, 14.4f, 15.0f, 16.0f, 18.0f, 19.2f, 20.0f}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.2f, 4.4f, 4.6f, 4.8f, 5.0f, 5.2f, 5.4f, 5.6f, 5.8f, 6.0f, 6.2f, 6.4f, 6.6f, 6.8f, 7.0f, 7.2f, 7.4f, 7.6f, 7.8f, 8.0f, 8.2f, 8.4f, 8.6f, 8.8f, 9.0f, 9.2f, 9.4f, 9.6f, 9.8f, 10.0f, 10.2f, 10.4f, 10.6f, 10.8f, 11.0f, 11.2f, 11.4f, 11.6f, 11.8f, 12.0f, 12.2f, 12.4f, 12.6f, 12.8f, 13.0f, 13.2f, 13.4f, 13.6f, 13.8f, 14.0f, 14.2f, 14.4f, 14.6f, 14.8f, 15.0f, 15.2f, 15.4f, 15.6f, 15.8f, 16.0f, 16.2f, 16.4f, 16.6f, 16.8f, 17.0f, 17.2f, 17.4f, 17.6f, 17.8f, 18.0f, 18.2f, 18.4f, 18.6f, 18.8f, 19.0f, 19.2f, 19.4f, 19.6f, 19.8f, 20.0f}, "Binning of the pT axis"};
  ConfigurableAxis binsEta{"binsEta", {150, -1.5f, 1.5f}, ""};
  ConfigurableAxis binsMass{"binsMass", {70, 1.3f, 2.0f}, "Invariant Mass (GeV/#it{c}^2)"};
  ConfigurableAxis binsMult{"binsMult", {105, 0.0f, 105.0f}, "mult_{FT0M}"};
  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2f, 0.2f}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2f, 0.2f}, ""};
  ConfigurableAxis binsTPCXrows{"binsTPCXrows", {100, 60, 160}, ""};
  ConfigurableAxis binsnSigma{"binsnSigma", {130, -6.5f, 6.5f}, ""};
  ConfigurableAxis binsnTPCSignal{"binsnTPCSignal", {1000, 0, 1000}, ""};
  ConfigurableAxis binsEtaPhi{"binsEtaPhi", {350, -3.5f, 3.5f}, ""};

  void init(framework::InitContext&)
  {
    colCuts.setCuts(configEvents.cfgEvtZvtx, /* configEvents.cfgEvtTriggerCheck */ false, configEvents.cfgEvtSel8, /*checkRun3*/ true, /*triggerTVXsel*/ false, configEvents.cfgEvtOccupancyInTimeRangeMax, configEvents.cfgEvtOccupancyInTimeRangeMin);

    colCuts.init(&histos);
    colCuts.setTriggerTVX(configEvents.cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(configEvents.cfgEvtNoTFBorderCut);
    colCuts.setApplyITSTPCvertex(configEvents.cfgEvtIsVertexITSTPC);
    colCuts.setApplyZvertexTimedifference(configEvents.cfgEvtIsGoodZvtxFT0vsPV);
    colCuts.setApplyPileupRejection(configEvents.cfgEvtNoSameBunchPileup);
    colCuts.setApplyNoITSROBorderCut(configEvents.cfgEvtNoITSROFrameBorderCut);
    colCuts.setApplyCollInTimeRangeStandard(configEvents.cfgEvtNoCollInTimeRangeStandard);
    colCuts.setApplyVertexTOFmatched(configEvents.cfgEvtIsVertexTOFmatched);
    colCuts.printCuts();

    // axes
    AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtQA{binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{binsEta, "#eta"};
    AxisSpec axisRap{binsEta, "#it{y}"};
    AxisSpec axisMassLambda1520{binsMass, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec axisMult{binsMult, "mult_{FT0M}"};
    AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};
    AxisSpec axisTPCXrow{binsTPCXrows, "#Xrows_{TPC}"};
    AxisSpec axisPIDQA{binsnSigma, "#sigma"};
    AxisSpec axisTPCSignal{binsnTPCSignal, ""};
    AxisSpec axisMClabel{6, -1.5f, 6.5f, "MC Label"};
    AxisSpec axisEtaPhi{binsEtaPhi, ""};
    AxisSpec axisPhi{350, 0, 7, "#Phi"};
    AxisSpec axisMultMix{configBkg.cfgMultPercentileBins, "Multiplicity Percentile"};
    AxisSpec axisVtxMix{configBkg.cfgVtxBins, "Vertex Z (cm)"};
    AxisSpec idxMCAxis = {26, -0.5f, 25.5f, "Index"};

    if (cFilladditionalQAeventPlots) {
      // event histograms
      if (doprocessData) {
        histos.add("QAevent/hEvents", "INEL>0 Events", HistType::kTH1F, {{2, 0.5f, 2.5f}});
        histos.add("QAevent/hPairsCounterSameE", "total valid no. of pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
        histos.add("QAevent/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0, 1000.0}});
      }
      // Test on Mixed event
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
      if (doprocessMCRec) {
        histos.add("QAevent/hEventsMC", "INEL>0 Events MC", HistType::kTH1F, {{2, 0.5f, 2.5f}});
      }
    }

    if (doprocessData) {
      // Track QA before cuts
      //  --- Track
      if (cFillTrackQA) {
        histos.add("QA/QAbefore/Track/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;{#sigma_{TOF}^{Kaon}};{#sigma_{TPC}^{Kaon}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPConly_Nsigma_ka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTH2F, {axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;{#sigma_{TOF}^{Proton}};{#sigma_{TPC}^{Proton}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Proton}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Proton}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/TPConly_Nsigma_pr", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Proton}};", {HistType::kTH2F, {axisPt, axisPIDQA}});
        histos.add("QA/QAbefore/Track/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
        histos.add("QA/QAbefore/Track/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
        histos.add("QA/QAbefore/Track/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
        histos.add("QA/QAbefore/Track/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
        histos.add("QA/QAbefore/Track/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      }
      if (cFillMultQA) {
        // Multiplicity correlation calibrations
        histos.add("MultCalib/centGloPVpr", "Centrality vs Global-Tracks", kTHnSparseF, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}, {500, 0, 5000, "PV tracks"}});
        histos.add("MultCalib/centGloPVka", "Centrality vs Global-Tracks", kTHnSparseF, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}, {500, 0, 5000, "PV tracks"}});
      }

      // PID QA after cuts
      //  --- Kaon
      histos.add("QA/QAafter/Kaon/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;{#sigma_{TOF}^{Kaon}};{#sigma_{TPC}^{Kaon}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTH2F, {axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Kaon/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Kaon/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Kaon/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Kaon/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Kaon/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Kaon/TPC_Signal_ka_all", "TPC Signal for Kaon;#it{p} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});
      histos.add("QA/QAafter/Kaon/TPCnclusterPhika", "TPC ncluster vs phi", kTHnSparseF, {{160, 0, 160, "TPC nCluster"}, {63, 0.0f, 6.28f, "#phi"}});

      //  --- Proton
      histos.add("QA/QAafter/Proton/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;{#sigma_{TOF}^{Proton}};{#sigma_{TPC}^{Proton}}", {HistType::kTH2F, {axisPIDQA, axisPIDQA}});
      histos.add("QA/QAafter/Proton/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Proton}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Proton}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Proton}};", {HistType::kTH2F, {axisPt, axisPIDQA}});
      histos.add("QA/QAafter/Proton/dcaZ", "DCA_{Z} distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm);", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Proton/dcaXY", "DCA_{XY} momentum distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Proton/TPC_CR", "# TPC Xrows distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Proton/pT", "pT distribution of Protons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Proton/eta", "#eta distribution of Protons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Proton/TPC_Signal_pr_all", "TPC Signal for Proton;#it{p} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});
      histos.add("QA/QAafter/Proton/TPCnclusterPhipr", "TPC ncluster vs phi", kTHnSparseF, {{160, 0, 160, "TPC nCluster"}, {63, 0.0f, 6.28f, "#phi"}});

      //  Mass QA 1D for quick check
      if (cFill1DQAs) {
        histos.add("Result/Data/lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/antilambda1520invmass", "Invariant mass of #Lambda(1520) K^{#mp}p^{#pm}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});   // K+ + Pr
        histos.add("Result/Data/lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}}); // K- + anti-Pr
      }
      // eta phi QA
      if (cFilldeltaEtaPhiPlots) {
        histos.add("QAbefore/deltaEta", "deltaEta of kaon and proton candidates", HistType::kTH1F, {axisEtaPhi});
        histos.add("QAbefore/deltaPhi", "deltaPhi of kaon and proton candidates", HistType::kTH1F, {axisEtaPhi});

        histos.add("QAafter/deltaEta", "deltaEta of kaon and proton candidates", HistType::kTH1F, {axisEtaPhi});
        histos.add("QAafter/deltaPhi", "deltaPhi of kaon and proton candidates", HistType::kTH1F, {axisEtaPhi});

        histos.add("QAafter/PhiPrafter", "Phi of  proton candidates", HistType::kTH1F, {axisPhi});
        histos.add("QAafter/PhiKaafter", "Phi of kaon  candidates", HistType::kTH1F, {axisPhi});
      }

      // 3d histogram
      histos.add("Result/Data/h3lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3antilambda1520invmass", "Invariant mass of #Lambda(1520) K^{#mp}p^{#pm}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});   // K+ + Pr
      histos.add("Result/Data/h3lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520}); // K- + anti-Pr
    }

    if (doprocessRotational) {
      if (cFill1DQAs) {
        histos.add("Result/Data/lambda1520InvMassRotation", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});       // K+ + Pr
        histos.add("Result/Data/antilambda1520InvMassRotation", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}}); // K- + anti-Pr
      }
      histos.add("Result/Data/h3lambda1520InvMassRotation", "Invariant mass of #Lambda(1520) rotation", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3antilambda1520InvMassRotation", "Invariant mass of #Lambda(1520) rotation", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
    }
    // Mixed event histograms
    if (doprocessME) {
      if (cFill1DQAs) {
        histos.add("Result/Data/lambda1520invmassME_UnlikeSign", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/antilambda1520invmassME_UnlikeSign", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
      }
      histos.add("Result/Data/h3lambda1520invmassME_UnlikeSign", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3antilambda1520invmassME_UnlikeSign", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});

      if (cFilladditionalMEPlots) {
        if (cFill1DQAs) {
          histos.add("Result/Data/lambda1520invmassME_LSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});   // K+ + Pr
          histos.add("Result/Data/lambda1520invmassME_LSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}}); // K- + anti-Pr
        }
        histos.add("Result/Data/h3lambda1520invmassME_LSPP", "Invariant mass of #Lambda(1520) mixed event Like Sign Method K^{#plus}p^{#plus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});   // K+ + Pr
        histos.add("Result/Data/h3lambda1520invmassME_LSMM", "Invariant mass of #Lambda(1520) mixed event Like Sign Method K^{#minus}p^{#minus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520}); // K- + anti-Pr
      }
    }

    // MC QA
    histos.add("Event/hMCEventIndices", "hMCEventIndices", kTH2D, {axisMult, idxMCAxis});
    if (doprocessMCGen) {
      histos.add("QA/MC/h2GenEtaPt_beforeanycut", " #eta-#it{p}_{T} distribution of Generated #Lambda(1520); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_beforeanycut", " #phi-y distribution of Generated #Lambda(1520); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});
      histos.add("QA/MC/h2GenEtaPt_afterEtaRapCut", " #eta-#it{p}_{T} distribution of Generated #Lambda(1520); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_afterEtaRapCut", " #phi-y distribution of Generated #Lambda(1520); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});
      histos.add("QA/MC/h2GenEtaPt_afterRapcut", " #phi-#it{p}_{T} distribution of Generated #Lambda(1520); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_afterRapcut", " #phi-y distribution of Generated #Lambda(1520); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});

      histos.add("Result/MC/Genlambda1520pt", "pT distribution of True MC #Lambda(1520)0", kTH3F, {axisMClabel, axisPt, axisMult});
      histos.add("Result/MC/Genantilambda1520pt", "pT distribution of True MC Anti-#Lambda(1520)0", kTH3F, {axisMClabel, axisPt, axisMult});
    }
    if (doprocessMCRec) {
      histos.add("QA/MC/h2RecoEtaPt_after", " #eta-#it{p}_{T} distribution of Reconstructed #Lambda(1520); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPt});
      histos.add("QA/MC/h2RecoPhiRapidity_after", " #phi-y distribution of Reconstructed #Lambda(1520); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});

      histos.add("QA/MC/trkDCAxy_pr", "DCAxy distribution of proton track candidates", HistType::kTHnSparseF, {axisPt, axisDCAxy});
      histos.add("QA/MC/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTHnSparseF, {axisPt, axisDCAxy});
      histos.add("QA/MC/trkDCAz_pr", "DCAz distribution of proton track candidates", HistType::kTHnSparseF, {axisPt, axisDCAz});
      histos.add("QA/MC/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTHnSparseF, {axisPt, axisDCAz});
      histos.add("QA/MC/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Proton}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/MC/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Proton}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/MC/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TOF}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});
      histos.add("QA/MC/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});{#sigma_{TPC}^{Kaon}};", {HistType::kTHnSparseF, {axisMult, axisPt, axisPIDQA}});

      histos.add("Result/MC/h3lambda1520Recoinvmass", "Invariant mass of Reconstructed MC #Lambda(1520)0", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/MC/h3antilambda1520Recoinvmass", "Invariant mass of Reconstructed MC Anti-#Lambda(1520)0", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
    }

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in Lambda1520analysisinpp:";
    histos.print();
  }

  float massKa = MassKaonCharged;
  float massPr = MassProton;

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
  bool ptDependentPidProton(const T& candidate)
  {
    auto vProtonTPCPIDpTintv = configPID.protonTPCPIDpTintv.value;
    vProtonTPCPIDpTintv.insert(vProtonTPCPIDpTintv.begin(), configTracks.cMinPtcut);
    auto vProtonTPCPIDcuts = configPID.protonTPCPIDcuts.value;
    auto vProtonTOFPIDpTintv = configPID.protonTOFPIDpTintv.value;
    auto vProtonTPCTOFCombinedpTintv = configPID.protonTPCTOFCombinedpTintv.value;
    auto vProtonTPCTOFCombinedPIDcuts = configPID.protonTPCTOFCombinedPIDcuts.value;
    auto vProtonTOFPIDcuts = configPID.protonTOFPIDcuts.value;

    float pt = candidate.pt();
    float ptSwitchToTOF = vProtonTPCPIDpTintv.back();
    float tpcNsigmaPr = candidate.tpcNSigmaPr();
    float tofNsigmaPr = candidate.tofNSigmaPr();

    bool tpcPIDPassed = false;

    // TPC PID (interval check)
    for (size_t i = 0; i < vProtonTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vProtonTPCPIDpTintv[i] && pt < vProtonTPCPIDpTintv[i + 1]) {
        if (std::abs(tpcNsigmaPr) < vProtonTPCPIDcuts[i])
          tpcPIDPassed = true;
      }
    }

    // TOF bypass option (for QA or MC)
    if (configPID.cByPassTOF) {
      return std::abs(tpcNsigmaPr) < vProtonTPCPIDcuts.back();
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
        for (size_t i = 0; i < vProtonTOFPIDpTintv.size(); ++i) {
          if (pt < vProtonTOFPIDpTintv[i]) {
            if (std::abs(tofNsigmaPr) < vProtonTOFPIDcuts[i] &&
                std::abs(tpcNsigmaPr) < vProtonTPCPIDcuts.back())
              return true;
          }
        }
      } else if (configPID.cPIDcutType == CircularType) {
        // Circular cut
        for (size_t i = 0; i < vProtonTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vProtonTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 =
              tpcNsigmaPr * tpcNsigmaPr +
              tofNsigmaPr * tofNsigmaPr;
            if (combinedSigma2 < vProtonTPCTOFCombinedPIDcuts[i] * vProtonTPCTOFCombinedPIDcuts[i])
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

  auto static constexpr MinPtforPionRejection = 1.0f;
  auto static constexpr MaxPtforPionRejection = 2.0f;
  auto static constexpr MaxnSigmaforPionRejection = 2.0f;

  template <typename T>
  bool rejectPion(const T& candidate)
  {
    if (candidate.pt() > MinPtforPionRejection && candidate.pt() < MaxPtforPionRejection && !candidate.hasTOF() && candidate.tpcNSigmaPi() < MaxnSigmaforPionRejection) {
      return false;
    }
    return true;
  }

  auto static constexpr MaxNoLambda1520Daughters = 2;

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
      // Trk1: Proton
      auto isTrk1hasTOF = trk1.hasTOF();
      auto trk1ptPr = trk1.pt();
      auto trk1etaPr = trk1.eta();
      auto trk1phiPr = trk1.phi();
      auto trk1NSigmaPrTPC = trk1.tpcNSigmaPr();
      auto trk1NSigmaPrTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPr() : -999.0f;

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
        deltaEta = std::abs(trk1etaPr - trk2etaKa);
        deltaPhi = std::abs(trk1phiPr - trk2phiKa);
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
          histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_pr_all"), centrality, trk1ptPr, trk1NSigmaPrTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_pr_all"), centrality, trk1ptPr, trk1NSigmaPrTOF);
            histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
          }
          if (!isTrk1hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_pr"), trk1ptPr, trk1NSigmaPrTPC);
          }
          histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTOF);
            histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
          }
          if (!isTrk2hasTOF) {
            histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_ka"), trk2ptKa, trk2NSigmaKaTPC);
          }

          histos.fill(HIST("QA/QAbefore/Track/dcaZ"), trk1ptPr, trk1.dcaZ());
          histos.fill(HIST("QA/QAbefore/Track/dcaXY"), trk1ptPr, trk1.dcaXY());
          histos.fill(HIST("QA/QAbefore/Track/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
          histos.fill(HIST("QA/QAbefore/Track/pT"), trk1ptPr);
          histos.fill(HIST("QA/QAbefore/Track/eta"), trk1etaPr);
        }
        if (cFilldeltaEtaPhiPlots) {
          histos.fill(HIST("QAbefore/deltaEta"), deltaEta);
          histos.fill(HIST("QAbefore/deltaPhi"), deltaPhi);
        }
      }

      //// Apply the pid selection
      if (crejectPion && rejectPion(trk2)) // to remove pion contamination from the kaon track
        continue;

      if (!ptDependentPidProton(trk1) || !ptDependentPidKaon(trk2))
        continue;

      //// QA plots after the selection
      if constexpr (IsData) { //  --- PID QA Proton
        histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_all"), centrality, trk1ptPr, trk1NSigmaPrTPC);
        histos.fill(HIST("QA/QAafter/Proton/TPC_Signal_pr_all"), trk1.tpcInnerParam(), trk1.tpcSignal());
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TOF_Nsigma_pr_all"), centrality, trk1ptPr, trk1NSigmaPrTOF);
          histos.fill(HIST("QA/QAafter/Proton/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly"), trk1ptPr, trk1NSigmaPrTPC);
        }
        histos.fill(HIST("QA/QAafter/Proton/dcaZ"), trk1ptPr, trk1.dcaZ());
        histos.fill(HIST("QA/QAafter/Proton/dcaXY"), trk1ptPr, trk1.dcaXY());
        histos.fill(HIST("QA/QAafter/Proton/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Proton/pT"), trk1ptPr);
        histos.fill(HIST("QA/QAafter/Proton/eta"), trk1etaPr);
        histos.fill(HIST("QA/QAafter/Proton/TPCnclusterPhipr"), trk1.tpcNClsFound(), trk1phiPr);

        //  --- PID QA Kaon
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTPC);
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Signal_ka_all"), trk2.tpcInnerParam(), trk2.tpcSignal());
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TOF_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAafter/Kaon/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly"), trk2ptKa, trk2NSigmaKaTPC);
        }
        histos.fill(HIST("QA/QAafter/Kaon/dcaZ"), trk2ptKa, trk2.dcaZ());
        histos.fill(HIST("QA/QAafter/Kaon/dcaXY"), trk2ptKa, trk2.dcaXY());
        histos.fill(HIST("QA/QAafter/Kaon/TPC_CR"), trk2ptKa, trk2.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Kaon/pT"), trk2ptKa);
        histos.fill(HIST("QA/QAafter/Kaon/eta"), trk2etaKa);
        histos.fill(HIST("QA/QAafter/Kaon/TPCnclusterPhika"), trk2.tpcNClsFound(), trk2phiKa);

        if (cFilldeltaEtaPhiPlots) {
          histos.fill(HIST("QAafter/PhiPrafter"), trk1phiPr);
          histos.fill(HIST("QAafter/PhiKaafter"), trk2phiKa);
          histos.fill(HIST("QAafter/deltaEta"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhi"), deltaPhi);
        }
      }

      //// Resonance reconstruction
      lDecayDaughter1 = LorentzVectorPtEtaPhiMass(trk1ptPr, trk1etaPr, trk1phiPr, massPr);
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
            if (configBkg.cfgRotPr) {
              ldaughterRot = LorentzVectorPtEtaPhiMass(trk1ptPr, trk1etaPr, trk1phiPr + theta, massPr);
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
                histos.fill(HIST("Result/Data/lambda1520InvMassRotation"), resonanceRotMass);
              }
              histos.fill(HIST("Result/Data/h3lambda1520InvMassRotation"), centrality, resonanceRotPt, resonanceRotMass);
            } else if (trk1.sign() > 0) {
              if (cFill1DQAs) {
                histos.fill(HIST("Result/Data/antilambda1520InvMassRotation"), resonanceRotMass);
              }
              histos.fill(HIST("Result/Data/h3antilambda1520InvMassRotation"), centrality, resonanceRotPt, resonanceRotMass);
            }
          }
        }
        if constexpr (IsData) {
          if (trk1.sign() < 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/lambda1520invmass"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmass"), centrality, resonancePt, resonanceMass);
          } else if (trk1.sign() > 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/antilambda1520invmass"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3antilambda1520invmass"), centrality, resonancePt, resonanceMass);
          }
        } else if (IsMix) {
          if (trk1.sign() < 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/lambda1520invmassME_UnlikeSign"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmassME_UnlikeSign"), centrality, resonancePt, resonanceMass);
          } else if (trk1.sign() > 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/antilambda1520invmassME_UnlikeSign"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3antilambda1520invmassME_UnlikeSign"), centrality, resonancePt, resonanceMass);
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
          while (motherstrk1.size() > MaxNoLambda1520Daughters) {
            motherstrk1.pop_back();
            mothersPDGtrk1.pop_back();
          }

          const auto& mctrk2 = trk2.mcParticle();
          if (mctrk2.has_mothers()) {
            motherstrk2 = getMothersIndeces(mctrk2);
            mothersPDGtrk2 = getMothersPDGCodes(mctrk2);
          }
          while (motherstrk2.size() > MaxNoLambda1520Daughters) {
            motherstrk2.pop_back();
            mothersPDGtrk2.pop_back();
          }

          if (std::abs(mctrk1.pdgCode()) != kProton || std::abs(mctrk2.pdgCode()) != kKPlus)
            continue;

          if (motherstrk1[0] != motherstrk2[0]) // Same mother
            continue;

          if (std::abs(mothersPDGtrk1[0]) != Pdg::kLambda1520_Py)
            continue;

          // LOGF(info, "mother trk1 id: %d, mother trk1: %d, trk1 id: %d, trk1 pdgcode: %d, mother trk2 id: %d, mother trk2: %d, trk2 id: %d, trk2 pdgcode: %d", motherstrk1[0], mothersPDGtrk1[0], trk1.globalIndex(), mctrk1.pdgCode(), motherstrk2[0], mothersPDGtrk2[0], trk2.globalIndex(), mctrk2.pdgCode());

          if (cUseEtacutMC && std::abs(lResonance.Eta()) > cEtacutMC) // eta cut
            continue;

          if (cUseRapcutMC && std::abs(resonanceRapidity) > configTracks.cfgCutRapidity) // rapidity cut
            continue;

          histos.fill(HIST("QA/MC/h2RecoEtaPt_after"), lResonance.Eta(), resonancePt);
          histos.fill(HIST("QA/MC/h2RecoPhiRapidity_after"), lResonance.Phi(), resonanceRapidity);

          // Track selection check.
          histos.fill(HIST("QA/MC/trkDCAxy_pr"), trk1ptPr, trk1.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAxy_ka"), trk2ptKa, trk2.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAz_pr"), trk1ptPr, trk1.dcaZ());
          histos.fill(HIST("QA/MC/trkDCAz_ka"), trk2ptKa, trk2.dcaZ());

          histos.fill(HIST("QA/MC/TPC_Nsigma_pr_all"), centrality, trk1ptPr, trk1NSigmaPrTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_pr_all"), centrality, trk1ptPr, trk1NSigmaPrTOF);
          }
          histos.fill(HIST("QA/MC/TPC_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_ka_all"), centrality, trk2ptKa, trk2NSigmaKaTOF);
          }

          // MC histograms
          if (mothersPDGtrk1[0] > 0) {
            histos.fill(HIST("Result/MC/h3lambda1520Recoinvmass"), centrality, resonancePt, resonanceMass);
          } else {
            histos.fill(HIST("Result/MC/h3antilambda1520Recoinvmass"), centrality, resonancePt, resonanceMass);
          }
        }
      } else {
        if constexpr (IsData) {
          // Like sign pair ++
          if (trk1.sign() > 0) {
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/lambda1520invmassLSPP"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSPP"), centrality, resonancePt, resonanceMass);
          } else { // Like sign pair --
            if (cFill1DQAs) {
              histos.fill(HIST("Result/Data/lambda1520invmassLSMM"), resonanceMass);
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSMM"), centrality, resonancePt, resonanceMass);
          }
        } else if (IsMix) {
          if (cFilladditionalMEPlots) {
            // Like sign pair ++
            if (trk1.sign() > 0) {
              if (cFill1DQAs) {
                histos.fill(HIST("Result/Data/lambda1520invmassME_LSPP"), resonanceMass);
              }
              histos.fill(HIST("Result/Data/h3lambda1520invmassME_LSPP"), centrality, resonancePt, resonanceMass);
            } else { // Like sign pair --
              if (cFill1DQAs) {
                histos.fill(HIST("Result/Data/lambda1520invmassME_LSMM"), resonanceMass);
              }
              histos.fill(HIST("Result/Data/h3lambda1520invmassME_LSMM"), centrality, resonancePt, resonanceMass);
            }
          }
        }
      }
    }
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks)
  {
    if (!colCuts.isSelected(collision)) // Default event selection
      return;

    if (cFilladditionalQAeventPlots) {
      histos.fill(HIST("QAevent/hEvents"), 1);
    }

    if (!collision.isInelGt0()) // <--
      return;

    if (cFilladditionalQAeventPlots) {
      histos.fill(HIST("QAevent/hEvents"), 2);
    }

    colCuts.fillQA(collision);

    fillHistograms<true, false, false, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Lambda1520analysisinpp, processData, "Process Event for data without partition", false);

  void processRotational(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    if (!colCuts.isSelected(collision, false)) // Default event selection
      return;

    if (!collision.isInelGt0()) // <--
      return;

    fillHistograms<false, true, false, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Lambda1520analysisinpp, processRotational, "Process Rotational Background", false);

  void processMCRec(MCEventCandidates::iterator const& collision,
                    aod::McCollisions const&,
                    MCTrackCandidates const& tracks, aod::McParticles const&)
  {
    if (!colCuts.isSelected(collision))
      return;

    if (cFilladditionalQAeventPlots) {
      histos.fill(HIST("QAevent/hEventsMC"), 1);
    }

    if (!collision.isInelGt0()) // <--
      return;

    if (cFilladditionalQAeventPlots) {
      histos.fill(HIST("QAevent/hEventsMC"), 2);
    }

    fillHistograms<false, false, true, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Lambda1520analysisinpp, processMCRec, "Process Event for MC Rec without partition", false);

  Partition<aod::McParticles> selectedMCParticles = (nabs(aod::mcparticle::pdgCode) == static_cast<int>(Pdg::kLambda1520_Py)); // Lambda(1520)

  void processMCGen(MCEventCandidates::iterator const& collision, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    bool isInAfterAllCuts = colCuts.isSelected(collision, false);
    bool inVtx10 = (std::abs(collision.mcCollision().posZ()) > configEvents.cfgEvtZvtx) ? false : true;
    bool isTriggerTVX = collision.selection_bit(aod::evsel::kIsTriggerTVX);
    bool isSel8 = collision.sel8();

    auto mcPartsAll = mcParticles.sliceBy(perMcCollision, collision.mcCollision().globalIndex());

    bool isTrueINELgt0 = pwglf::isINELgt0mc(mcPartsAll, pdg);
    // bool isTrueINELgt0 = collision.isInelGt0();

    auto centrality = centEst(collision);

    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());

    // Not related to the real collisions
    for (const auto& part : mcParts) { // loop over all MC particles

      std::vector<int> daughterPDGs;
      if (part.has_daughters()) {
        auto daughter01 = mcParticles.rawIteratorAt(part.daughtersIds()[0] - mcParticles.offset());
        auto daughter02 = mcParticles.rawIteratorAt(part.daughtersIds()[1] - mcParticles.offset());
        daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
      } else {
        daughterPDGs = {-1, -1};
      }

      bool pass1 = std::abs(daughterPDGs[0]) == kKPlus || std::abs(daughterPDGs[1]) == kKPlus;   // At least one decay to Kaon
      bool pass2 = std::abs(daughterPDGs[0]) == kProton || std::abs(daughterPDGs[1]) == kProton; // At least one decay to Proton

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

      if (cUseEtacutMC && std::abs(part.eta()) > cEtacutMC) // eta cut
        continue;

      histos.fill(HIST("QA/MC/h2GenEtaPt_afterEtaRapCut"), part.eta(), part.pt());
      histos.fill(HIST("QA/MC/h2GenPhiRapidity_afterEtaRapCut"), part.phi(), part.y());

      // without any event selection
      if (part.pdgCode() > 0)
        histos.fill(HIST("Result/MC/Genlambda1520pt"), 0, part.pt(), centrality);
      else
        histos.fill(HIST("Result/MC/Genantilambda1520pt"), 0, part.pt(), centrality);

      if (inVtx10) // INEL10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 1, part.pt(), centrality);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 1, part.pt(), centrality);
      }
      if (inVtx10 && isSel8) // INEL>10, vtx10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 2, part.pt(), centrality);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 2, part.pt(), centrality);
      }
      if (inVtx10 && isTriggerTVX) // vtx10, TriggerTVX
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 3, part.pt(), centrality);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 3, part.pt(), centrality);
      }
      if (isInAfterAllCuts) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 4, part.pt(), centrality);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 4, part.pt(), centrality);
      }
      if (isInAfterAllCuts && isTrueINELgt0) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 5, part.pt(), centrality);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 5, part.pt(), centrality);
      }
    }

    // QA for Trigger efficiency
    histos.fill(HIST("Event/hMCEventIndices"), centrality, Inel);
    if (inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Inel10);
    if (isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Inelg0);
    if (inVtx10 && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Inelg010);

    // TVX MB trigger
    if (isTriggerTVX)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Trig);
    if (isTriggerTVX && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Trig10);
    if (isTriggerTVX && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, TrigINELg0);
    if (isTriggerTVX && isTrueINELgt0 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, TrigINELg010);

    // Sel8 event selection
    if (isSel8)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Sel8);
    if (isSel8 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Sel810);
    if (isSel8 && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Sel8INELg0);
    if (isSel8 && isTrueINELgt0 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, Sel8INELg010);

    // CollisionCuts selection
    if (isInAfterAllCuts)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, AllCuts);
    if (isInAfterAllCuts && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, AllCuts10);
    if (isInAfterAllCuts && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, AllCutsINELg0);
    if (isInAfterAllCuts && isTrueINELgt0 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, AllCutsINELg010);
  }
  PROCESS_SWITCH(Lambda1520analysisinpp, processMCGen, "Process Event for MC only", false);

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

      if (!colCuts.isSelected(collision1, false)) // Default event selection
        return;

      if (!colCuts.isSelected(collision2, false)) // Default event selection
        return;

      if (!collision1.isInelGt0()) // <--
        return;

      if (!collision2.isInelGt0()) // <--
        return;

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
  PROCESS_SWITCH(Lambda1520analysisinpp, processME, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Lambda1520analysisinpp>(cfgc)};
}
