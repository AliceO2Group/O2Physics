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

/// \file lstaranalysis.cxx
/// \brief This standalone task reconstructs track-track decay of lambda(1520) resonance candidate
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

// 1. Own header (doesn't exist)

// 2. C system headers (none)

// 3. C++ system headers
#include <string>

// 4. Other includes: O2 framework, ROOT, etc.
#include "PWGLF/Utils/collisionCuts.h"

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
#include "TVector3.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using LorentzVectorPtEtaPhiMass = ROOT::Math::PtEtaPhiMVector;

enum {
  kINEL = 1,
  kINEL10,
  kINELg0,
  kINELg010,
  kTrig,
  kTrig10,
  kTrigINELg0,
  kTrigINELg010,
  kSel8,
  kSel810,
  kSel8INELg0,
  kSel8INELg010,
  kAllCuts,
  kAllCuts10,
  kAllCutsINELg0,
  kAllCutsINELg010,
};

struct Lstaranalysis {
  // Define slice per Resocollision
  SliceCache cache;
  Preslice<Tracks> perCollision = o2::aod::track::collisionId;
  Preslice<McParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<framework::O2DatabasePDG> pdg;

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;

  Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<int> cfgEvtOccupancyInTimeRangeMax{"cfgEvtOccupancyInTimeRangeMax", -1, "Evt sel: maximum track occupancy"};
  Configurable<int> cfgEvtOccupancyInTimeRangeMin{"cfgEvtOccupancyInTimeRangeMin", -1, "Evt sel: minimum track occupancy"};
  Configurable<bool> cfgEvtTriggerCheck{"cfgEvtTriggerCheck", false, "Evt sel: check for trigger"};
  Configurable<bool> cfgEvtOfflineCheck{"cfgEvtOfflineCheck", true, "Evt sel: check for offline selection"};
  Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", false, "Evt sel: triggerTVX selection (MB)"};
  Configurable<bool> cfgEvtTFBorderCut{"cfgEvtTFBorderCut", false, "Evt sel: apply TF border cut"};
  Configurable<bool> cfgEvtUseITSTPCvertex{"cfgEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
  Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", false, "Evt sel: apply Z-vertex time difference"};
  Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", false, "Evt sel: apply pileup rejection"};
  Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
  Configurable<bool> cfgEvtCollInTimeRangeStandard{"cfgEvtCollInTimeRangeStandard", false, "Evt sel: apply NoCollInTimeRangeStandard"};

  Configurable<float> cfgEventCentralityMin{"cfgEventCentralityMin", 0.0f, "Event sel: minimum centrality"};
  Configurable<float> cfgEventCentralityMax{"cfgEventCentralityMax", 100.0f, "Event sel: maximum centrality"};

  // Configurables
  // Pre-selection Track cuts
  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15f, "Minimal pT for tracks"};
  Configurable<float> cMinTPCNClsFound{"cMinTPCNClsFound", 120, "minimum TPCNClsFound value for good track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
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

  /// PID Selections
  Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};                       // By pass TOF PID selection
  Configurable<int> cPIDcutType{"cPIDcutType", 2, "cPIDcutType = 1 for square cut, 2 for circular cut"}; // By pass TOF PID selection

  // Kaon
  Configurable<std::vector<float>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {0.5}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<float>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<float>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<float>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
  Configurable<std::vector<float>> kaonTPCTOFCombinedpTintv{"kaonTPCTOFCombinedpTintv", {999.}, "pT intervals for Kaon TPC-TOF PID cuts"};
  Configurable<std::vector<float>> kaonTPCTOFCombinedPIDcuts{"kaonTPCTOFCombinedPIDcuts", {2}, "nSigma list for Kaon TPC-TOF PID cuts"};

  // Proton
  Configurable<std::vector<float>> protonTPCPIDpTintv{"protonTPCPIDpTintv", {0.9}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<float>> protonTPCPIDcuts{"protonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<float>> protonTOFPIDpTintv{"protonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<float>> protonTOFPIDcuts{"protonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
  Configurable<std::vector<float>> protonTPCTOFCombinedpTintv{"protonTPCTOFCombinedpTintv", {999.}, "pT intervals for Proton TPC-TOF PID cuts"};
  Configurable<std::vector<float>> protonTPCTOFCombinedPIDcuts{"protonTPCTOFCombinedPIDcuts", {2}, "nSigma list for Proton TPC-TOF PID cuts"};

  // Additional purity check
  Configurable<bool> crejectPion{"crejectPion", false, "Switch to turn on/off pion contamination"};
  Configurable<bool> cApplyOpeningAngle{"cApplyOpeningAngle", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<float> cMinOpeningAngle{"cMinOpeningAngle", 1.4, "Maximum deltaEta between daughters"};
  Configurable<float> cMaxOpeningAngle{"cMaxOpeningAngle", 2.4, "Maximum deltaPhi between daughters"};

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - multiplicity"};

  // Rotational background
  Configurable<bool> isCalcRotBkg{"isCalcRotBkg", true, "Calculate rotational background"};
  Configurable<int> rotationalcut{"rotationalcut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  Configurable<int> cNofRotations{"cNofRotations", 3, "Number of random rotations in the rotational background"};

  // MC selection cut
  Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};
  Configurable<float> cEtacutMC{"cEtacutMC", 0.5, "MC eta cut"};
  Configurable<bool> cUseRapcutMC{"cUseRapcutMC", true, "MC eta cut"};
  Configurable<bool> cUseEtacutMC{"cUseEtacutMC", true, "MC eta cut"};

  // cuts on mother
  Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", false, "Enable additional cuts on mother"};
  Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 10.0, "Maximum pt of mother cut"};
  Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 3.0, "Maximum Minv of mother cut"};
  Configurable<float> cMaxDeltaEtaCut{"cMaxDeltaEtaCut", 0.7, "Maximum deltaEta between daughters"};
  Configurable<float> cMaxDeltaPhiCut{"cMaxDeltaPhiCut", 1.5, "Maximum deltaPhi between daughters"};

  // switches
  Configurable<bool> cFillMultQA{"cFillMultQA", false, "Turn on/off additional QA plots"};
  Configurable<bool> cFilladditionalQAeventPlots{"cFilladditionalQAeventPlots", false, "Additional QA event plots"};
  Configurable<bool> cFilladditionalMEPlots{"cFilladditionalMEPlots", false, "Additional Mixed event plots"};
  Configurable<bool> cFilldeltaEtaPhiPlots{"cFilldeltaEtaPhiPlots", false, "Enamble additional cuts on daughters"};
  Configurable<bool> cFillinvmass1DPlots{"cFillinvmass1DPlots", false, "Invariant mass 1D"};
  Configurable<int> multEstimator{"multEstimator", 0, "Select multiplicity estimator: 0 - FT0M, 1 - FT0A, 2 - FT0C"};

  Configurable<int> cfgCentEst{"cfgCentEst", 2, "Centrality estimator, 1: FT0C, 2: FT0M"};

  TRandom* rn = new TRandom();

  // Pre-filters for efficient process
  // Filter tofPIDFilter = aod::track::tofExpMom < 0.f || ((aod::track::tofExpMom > 0.f) && ((nabs(aod::pidtof::tofNSigmaPi) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaKa) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaPr) < pidnSigmaPreSelectionCut))); // TOF
  // Filter tpcPIDFilter = nabs(aod::pidtpc::tpcNSigmaPi) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaKa) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaPr) < pidnSigmaPreSelectionCut; // TPC
  /* Filter trackFilter = (trackSelection == 0) ||
                       ((trackSelection == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks)); */
  Filter trackEtaFilter = nabs(aod::track::eta) < cfgCutEta; // Eta cut

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>>;

  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Filtered<soa::Join<TrackCandidates, aod::McTrackLabels>>;

  /// Figures
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 4.6, 4.8, 4.9, 5.0, 5.5, 5.6, 6.0, 6.4, 6.5, 7.0, 7.2, 8.0, 9.0, 9.5, 9.6, 10.0, 11.0, 11.5, 12.0, 13.0, 14.0, 14.4, 15.0, 16.0, 18.0, 19.2, 20.}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 12.2, 12.4, 12.6, 12.8, 13, 13.2, 13.4, 13.6, 13.8, 14, 14.2, 14.4, 14.6, 14.8, 15, 15.2, 15.4, 15.6, 15.8, 16, 16.2, 16.4, 16.6, 16.8, 17, 17.2, 17.4, 17.6, 17.8, 18, 18.2, 18.4, 18.6, 18.8, 19, 19.2, 19.4, 19.6, 19.8, 20}, "Binning of the pT axis"};
  ConfigurableAxis binsEta{"binsEta", {150, -1.5, 1.5}, ""};
  ConfigurableAxis binsMass{"binsMass", {70, 1.3, 2.0}, "Invariant Mass (GeV/#it{c}^2)"};
  ConfigurableAxis binsMult{"binsMult", {105, 0.0, 105.0}, "mult_{FT0M}"};
  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsTPCXrows{"binsTPCXrows", {100, 60, 160}, ""};
  ConfigurableAxis binsnSigma{"binsnSigma", {130, -6.5, 6.5}, ""};
  ConfigurableAxis binsnTPCSignal{"binsnTPCSignal", {1000, 0, 1000}, ""};
  ConfigurableAxis binsEtaPhi{"binsEtaPhi", {350, -3.5, 3.5}, ""};

  float centrality;

  void init(framework::InitContext&)
  {
    centrality = -999;

    colCuts.setCuts(cfgEvtZvtx, cfgEvtTriggerCheck, cfgEvtOfflineCheck, /*checkRun3*/ true, /*triggerTVXsel*/ false, cfgEvtOccupancyInTimeRangeMax, cfgEvtOccupancyInTimeRangeMin);

    colCuts.init(&histos);
    colCuts.setTriggerTVX(cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(cfgEvtTFBorderCut);
    colCuts.setApplyITSTPCvertex(cfgEvtUseITSTPCvertex);
    colCuts.setApplyZvertexTimedifference(cfgEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(cfgEvtPileupRejection);
    colCuts.setApplyNoITSROBorderCut(cfgEvtNoITSROBorderCut);
    colCuts.setApplyCollInTimeRangeStandard(cfgEvtCollInTimeRangeStandard);
    colCuts.printCuts();

    // axes
    AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtQA{binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{binsEta, "#eta"};
    AxisSpec axisRap{binsEta, "#it{y}"};
    AxisSpec axisMassLambda1520{binsMass, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec axisMult{binsMult, "mult_{V0M}"};
    AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};
    AxisSpec axisTPCXrow{binsTPCXrows, "#Xrows_{TPC}"};
    AxisSpec axisPIDQA{binsnSigma, "#sigma"};
    AxisSpec axisTPCSignal{binsnTPCSignal, ""};
    AxisSpec axisMClabel{6, -1.5, 5.5, "MC Label"};
    AxisSpec axisEtaPhi{binsEtaPhi, ""};
    AxisSpec axisPhi{350, 0, 7, "#Phi"};
    AxisSpec axisMultMix{cfgMultBins, "Multiplicity"};
    AxisSpec axisVtxMix{cfgVtxBins, "Vertex Z (cm)"};
    AxisSpec idxMCAxis = {26, -0.5, 25.5, "Index"};

    if (cFilladditionalQAeventPlots) {
      // event histograms
      if (doprocessData) {
        histos.add("TestME/hPairsCounterSameE", "tot n pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
        histos.add("QAevent/hEvtCounterSameE", "Number of analyzed Same Events", HistType::kTH1F, {{1, 0.5, 1.5}});
        histos.add("QAevent/hVertexZSameE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
        histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
        histos.add("TestME/hCollisionIndexSameE", "coll index sameE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
        histos.add("TestME/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      }
      // Test on Mixed event
      if (doprocessME) {

        // Histograms for Mixed Event Pool characteristics
        histos.add("QAevent/hMixPool_VtxZ", "Mixed Event Pool: Vertex Z;Vtx Z (cm);Counts", HistType::kTH1F, {axisVtxMix});
        histos.add("QAevent/hMixPool_Multiplicity", "Mixed Event Pool: Multiplicity;Multiplicity;Counts", HistType::kTH1F, {axisMultMix});
        histos.add("QAevent/hMixPool_VtxZ_vs_Multiplicity", "Mixed Event Pool: Vertex Z vs Multiplicity;Counts", HistType::kTH2F, {axisVtxMix, axisMultMix});

        histos.add("TestME/hPairsCounterMixedE", "tot n pairs mixedE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
        histos.add("QAevent/hEvtCounterMixedE", "Number of analyzed Mixed Events", HistType::kTH1F, {{1, 0.5, 1.5}});
        histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
        histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
        histos.add("TestME/hCollisionIndexMixedE", "coll index mixedE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
        histos.add("TestME/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      }
    }

    if (doprocessData) {
      // Track QA before cuts
      //  --- Track
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
      histos.add("QA/QAafter/Kaon/TPCnclusterPhika", "TPC ncluster vs phi", kTHnSparseF, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});

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
      histos.add("QA/QAafter/Proton/TPCnclusterPhipr", "TPC ncluster vs phi", kTHnSparseF, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});

      //  Mass QA 1D for quick check
      if (cFillinvmass1DPlots) {
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

        histos.add("QAafter/deltaEtaafter", "deltaEta of kaon and proton candidates", HistType::kTH1F, {axisEtaPhi});
        histos.add("QAafter/deltaPhiafter", "deltaPhi of kaon and proton candidates", HistType::kTH1F, {axisEtaPhi});
        histos.add("QAafter/EtaPrafter", "Eta of  proton candidates", HistType::kTH1F, {axisEta});
        histos.add("QAafter/PhiPrafter", "Phi of  proton candidates", HistType::kTH1F, {axisPhi});
        histos.add("QAafter/EtaKaafter", "Eta of kaon  candidates", HistType::kTH1F, {axisEta});
        histos.add("QAafter/PhiKaafter", "Phi of kaon  candidates", HistType::kTH1F, {axisPhi});
      }

      if (isCalcRotBkg) {
        histos.add("Result/Data/h3lambda1520InvMassRotation", "Invariant mass of #Lambda(1520) rotation", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      }

      // 3d histogram
      histos.add("Result/Data/h3lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3antilambda1520invmass", "Invariant mass of #Lambda(1520) K^{#mp}p^{#pm}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});   // K+ + Pr
      histos.add("Result/Data/h3lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520}); // K- + anti-Pr
    }
    if (doprocessME) {
      if (cFillinvmass1DPlots) {
        histos.add("Result/Data/lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
      }
      histos.add("Result/Data/h3lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});

      if (cFilladditionalMEPlots) {
        histos.add("Result/Data/h3lambda1520invmassME_DS", "Invariant mass of #Lambda(1520) mixed event DS", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassME_DSAnti", "Invariant mass of #Lambda(1520) mixed event DSAnti", kTHnSparseF, {axisMult, axisPt, axisMassLambda1520});
      }
    }

    // MC QA
    histos.add("Event/hMCEventIndices", "hMCEventIndices", kTH2D, {axisMult, idxMCAxis});
    if (doprocessMCTrue) {
      histos.add("QA/MC/h2GenEtaPt_beforeanycut", " #eta-#it{p}_{T} distribution of Generated #Lambda(1520); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_beforeanycut", " #phi-y distribution of Generated #Lambda(1520); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});
      histos.add("QA/MC/h2GenEtaPt_afterEtaRapCut", " #eta-#it{p}_{T} distribution of Generated #Lambda(1520); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_afterEtaRapCut", " #phi-y distribution of Generated #Lambda(1520); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});
      histos.add("QA/MC/h2GenEtaPt_afterRapcut", " #phi-#it{p}_{T} distribution of Generated #Lambda(1520); #eta;  #it{p}_{T}; Counts;", HistType::kTHnSparseF, {axisEta, axisPtQA});
      histos.add("QA/MC/h2GenPhiRapidity_afterRapcut", " #phi-y distribution of Generated #Lambda(1520); #phi; y; Counts;", HistType::kTHnSparseF, {axisPhi, axisRap});

      histos.add("Result/MC/Genlambda1520pt", "pT distribution of True MC #Lambda(1520)0", kTHnSparseF, {axisMClabel, axisPt, axisMult});
      histos.add("Result/MC/Genantilambda1520pt", "pT distribution of True MC Anti-#Lambda(1520)0", kTHnSparseF, {axisMClabel, axisPt, axisMult});
    }
    if (doprocessMC) {
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
    LOG(info) << "Size of the histograms in LstarAnalysis:";
    histos.print();
  }

  float massKa = MassKaonCharged;
  float massPr = MassProton;

  int kLambda1520PDG = static_cast<int>(102134); // PDG code for Lambda(1520)

  template <typename CollisionType>
  float getCentrality(CollisionType const& collision)
  {
    if (cfgCentEst == static_cast<int>(1)) {
      return collision.multFT0C();
    } else if (cfgCentEst == static_cast<int>(2)) {
      return collision.multFT0M();
    } else {
      return -999;
    }
  }

  // Centralicity estimator selection
  template <typename ResoColl>
  float centEst(ResoColl ResoEvents)
  {
    float returnValue = -999.0;
    switch (multEstimator) {
      case 0:
        returnValue = ResoEvents.centFT0M();
        break;
      case 1:
        returnValue = ResoEvents.centFT0A();
        break;
      case 2:
        returnValue = ResoEvents.centFT0C();
        break;
      default:
        returnValue = ResoEvents.centFT0M();
        break;
    }
    return returnValue;
  }

  // Check if the collision is INEL>0
  template <typename MCColl, typename MCPart>
  bool isTrueINEL0(MCColl const& /*mccoll*/, MCPart const& mcparts)
  {
    for (auto const& mcparticle : mcparts) {
      if (!mcparticle.isPhysicalPrimary())
        continue;
      auto p = pdg->GetParticle(mcparticle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          if (std::abs(mcparticle.eta()) < 1)
            return true;
        }
      }
    }
    return false;
  }

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (cDCAr7SigCut) {
      if (std::abs(track.dcaXY()) > (0.004f + 0.0130f / (track.pt()))) // 7 - Sigma cut
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
        return false;
    }
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (cTPCNClsFound && (track.tpcNClsFound() < cMinTPCNClsFound))
      return false;
    if (track.tpcNClsCrossedRows() < cfgMinCrossedRows)
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;

    return true;
  }

  // LOGF(info, "AFTER: pt: %f, hasTOF: %d, TPCSigma: %f, TOFSigma: %f, boolTPC: %d, boolTOF: %d, bool: %d", pt, candidate.hasTOF(),
  //    candidate.tpcNSigmaPr(), candidate.tofNSigmaPr(), tpcPIDPassed, tofPIDPassed, tpcPIDPassed || tofPIDPassed);

  template <typename T>
  bool pTdependentPIDProton(const T& candidate)
  {
    auto vProtonTPCPIDpTintv = static_cast<std::vector<float>>(protonTPCPIDpTintv);
    vProtonTPCPIDpTintv.insert(vProtonTPCPIDpTintv.begin(), cMinPtcut);
    auto vProtonTPCPIDcuts = static_cast<std::vector<float>>(protonTPCPIDcuts);
    auto vProtonTOFPIDpTintv = static_cast<std::vector<float>>(protonTOFPIDpTintv);
    auto vProtonTPCTOFCombinedpTintv = static_cast<std::vector<float>>(protonTPCTOFCombinedpTintv);
    auto vProtonTPCTOFCombinedPIDcuts = static_cast<std::vector<float>>(protonTPCTOFCombinedPIDcuts);
    auto vProtonTOFPIDcuts = static_cast<std::vector<float>>(protonTOFPIDcuts);

    float pt = candidate.pt();
    float ptSwitchToTOF = vProtonTPCPIDpTintv.back();

    bool tpcPIDPassed = false;

    // TPC PID (interval check)
    for (size_t i = 0; i < vProtonTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vProtonTPCPIDpTintv[i] && pt < vProtonTPCPIDpTintv[i + 1]) {
        if (std::abs(candidate.tpcNSigmaPr()) < vProtonTPCPIDcuts[i])
          tpcPIDPassed = true;
      }
    }

    // TOF bypass option (for QA or MC)
    if (cByPassTOF) {
      return std::abs(candidate.tpcNSigmaPr()) < vProtonTPCPIDcuts.back();
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
      if (cPIDcutType == 1) {
        // Rectangular cut
        for (size_t i = 0; i < vProtonTOFPIDpTintv.size(); ++i) {
          if (pt < vProtonTOFPIDpTintv[i]) {
            if (std::abs(candidate.tofNSigmaPr()) < vProtonTOFPIDcuts[i] &&
                std::abs(candidate.tpcNSigmaPr()) < vProtonTPCPIDcuts.back())
              return true;
          }
        }
      } else if (cPIDcutType == 2) {
        // Circular cut
        for (size_t i = 0; i < vProtonTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vProtonTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 =
              candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr() +
              candidate.tofNSigmaPr() * candidate.tofNSigmaPr();
            if (combinedSigma2 < vProtonTPCTOFCombinedPIDcuts[i] * vProtonTPCTOFCombinedPIDcuts[i])
              return true;
          }
        }
      }
    }

    return false;
  }

  template <typename T>
  bool pTdependentPIDKaon(const T& candidate)
  {
    auto vKaonTPCPIDpTintv = static_cast<std::vector<float>>(kaonTPCPIDpTintv);
    vKaonTPCPIDpTintv.insert(vKaonTPCPIDpTintv.begin(), cMinPtcut);
    auto vKaonTPCPIDcuts = static_cast<std::vector<float>>(kaonTPCPIDcuts);
    auto vKaonTOFPIDpTintv = static_cast<std::vector<float>>(kaonTOFPIDpTintv);
    auto vKaonTPCTOFCombinedpTintv = static_cast<std::vector<float>>(kaonTPCTOFCombinedpTintv);
    auto vKaonTPCTOFCombinedPIDcuts = static_cast<std::vector<float>>(kaonTPCTOFCombinedPIDcuts);
    auto vKaonTOFPIDcuts = static_cast<std::vector<float>>(kaonTOFPIDcuts);

    float pt = candidate.pt();
    float ptSwitchToTOF = vKaonTPCPIDpTintv.back();

    bool tpcPIDPassed = false;

    // TPC PID interval-based check
    for (size_t i = 0; i < vKaonTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vKaonTPCPIDpTintv[i] && pt < vKaonTPCPIDpTintv[i + 1]) {
        if (std::abs(candidate.tpcNSigmaKa()) < vKaonTPCPIDcuts[i]) {
          tpcPIDPassed = true;
          break;
        }
      }
    }

    // TOF bypass option
    if (cByPassTOF) {
      return std::abs(candidate.tpcNSigmaKa()) < vKaonTPCPIDcuts.back();
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
      if (cPIDcutType == 1) {
        // Rectangular cut
        for (size_t i = 0; i < vKaonTOFPIDpTintv.size(); ++i) {
          if (pt < vKaonTOFPIDpTintv[i]) {
            if (std::abs(candidate.tofNSigmaKa()) < vKaonTOFPIDcuts[i] &&
                std::abs(candidate.tpcNSigmaKa()) < vKaonTPCPIDcuts.back()) {
              return true;
            }
          }
        }
      } else if (cPIDcutType == 2) {
        // Circular cut
        for (size_t i = 0; i < vKaonTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vKaonTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 = candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() +
                                   candidate.tofNSigmaKa() * candidate.tofNSigmaKa();
            if (combinedSigma2 < vKaonTPCTOFCombinedPIDcuts[i] * vKaonTPCTOFCombinedPIDcuts[i]) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  template <typename T>
  bool rejectPion(const T& candidate)
  {
    if (candidate.pt() > static_cast<float>(1.0) && candidate.pt() < static_cast<float>(2.0) && !candidate.hasTOF() && candidate.tpcNSigmaPi() < static_cast<float>(2)) {
      return false;
    }
    return true;
  }

  template <bool IsData, bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.centFT0M();

    // Multiplicity correlation calibration plots
    if (cFillMultQA) {
      if constexpr (IsData) {
        histos.fill(HIST("MultCalib/centGloPVpr"), multiplicity, dTracks1.size(), collision.multNTracksPV());
        histos.fill(HIST("MultCalib/centGloPVka"), multiplicity, dTracks2.size(), collision.multNTracksPV());
      }
    }

    if (cFilladditionalQAeventPlots) {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), multiplicity);
        histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksSameE"), dTracks1.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), multiplicity);
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), dTracks1.size());
      }
    }
    // LOG(info) << "After pass, Collision index:" << collision.index() << "multiplicity: " << collision.centFT0M() << endl;

    LorentzVectorPtEtaPhiMass lDecayDaughter1, lDecayDaughter2, lResonance, ldaughterRot, lresonanceRot;

    for (const auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      if (cFilladditionalQAeventPlots) {
        if constexpr (IsData) {
          histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
        } else if (IsMix) {
          histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
        }
      }

      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      //// Initialize variables
      // Trk1: Proton, Trk2: Kaon
      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();

      auto trk1ptPr = trk1.pt();
      auto trk1NSigmaPrTPC = trk1.tpcNSigmaPr();
      auto trk1NSigmaPrTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPr() : -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

      auto deltaEta = std::abs(trk1.eta() - trk2.eta());
      auto deltaPhi = std::abs(trk1.phi() - trk2.phi());
      deltaPhi = (deltaPhi > constants::math::PI) ? (constants::math::TwoPI - deltaPhi) : deltaPhi;

      //// QA plots before the selection
      //  --- Track QA all
      if constexpr (IsData) {
        histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_pr"), trk1ptPr, trk1NSigmaPrTPC);
        }
        histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TPConly_Nsigma_ka"), trk2ptKa, trk2NSigmaKaTPC);
        }

        histos.fill(HIST("QA/QAbefore/Track/dcaZ"), trk1ptPr, trk1.dcaZ());
        histos.fill(HIST("QA/QAbefore/Track/dcaXY"), trk1ptPr, trk1.dcaXY());
        histos.fill(HIST("QA/QAbefore/Track/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAbefore/Track/pT"), trk1ptPr);
        histos.fill(HIST("QA/QAbefore/Track/eta"), trk1.eta());
        if (cFilldeltaEtaPhiPlots) {
          histos.fill(HIST("QAbefore/deltaEta"), deltaEta);
          histos.fill(HIST("QAbefore/deltaPhi"), deltaPhi);
        }
      }

      //// Apply the pid selection
      if (crejectPion && rejectPion(trk2))
        continue;

      if (!pTdependentPIDProton(trk1) || !pTdependentPIDKaon(trk2))
        continue;

      //// QA plots after the selection
      if constexpr (IsData) { //  --- PID QA Proton
        histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
        histos.fill(HIST("QA/QAafter/Proton/TPC_Signal_pr_all"), trk1.tpcInnerParam(), trk1.tpcSignal());
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          histos.fill(HIST("QA/QAafter/Proton/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly"), trk1ptPr, trk1NSigmaPrTPC);
        }
        histos.fill(HIST("QA/QAafter/Proton/dcaZ"), trk1ptPr, trk1.dcaZ());
        histos.fill(HIST("QA/QAafter/Proton/dcaXY"), trk1ptPr, trk1.dcaXY());
        histos.fill(HIST("QA/QAafter/Proton/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Proton/pT"), trk1ptPr);
        histos.fill(HIST("QA/QAafter/Proton/eta"), trk1.eta());
        histos.fill(HIST("QA/QAafter/Proton/TPCnclusterPhipr"), trk1.tpcNClsFound(), trk1.phi());

        //  --- PID QA Kaon
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Signal_ka_all"), trk2.tpcInnerParam(), trk2.tpcSignal());
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAafter/Kaon/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly"), trk2ptKa, trk2NSigmaKaTPC);
        }
        histos.fill(HIST("QA/QAafter/Kaon/dcaZ"), trk2ptKa, trk2.dcaZ());
        histos.fill(HIST("QA/QAafter/Kaon/dcaXY"), trk2ptKa, trk2.dcaXY());
        histos.fill(HIST("QA/QAafter/Kaon/TPC_CR"), trk2ptKa, trk2.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Kaon/pT"), trk2ptKa);
        histos.fill(HIST("QA/QAafter/Kaon/eta"), trk2.eta());
        histos.fill(HIST("QA/QAafter/Kaon/TPCnclusterPhika"), trk2.tpcNClsFound(), trk2.phi());

        if (cFilldeltaEtaPhiPlots) {
          histos.fill(HIST("QAafter/deltaEta"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhi"), deltaPhi);
        }
      }

      // Apply kinematic opening angle cut
      if (cApplyOpeningAngle) {
        TVector3 v1(trk1.px(), trk1.py(), trk1.pz());
        TVector3 v2(trk2.px(), trk2.py(), trk2.pz());
        float alpha = v1.Angle(v2);
        if (alpha > cMinOpeningAngle && alpha < cMaxOpeningAngle)
          continue;
      }

      //// Resonance reconstruction
      lDecayDaughter1 = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(), massPr);
      lDecayDaughter2 = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if constexpr (IsData || IsMix) {
        // Rapidity cut
        if (std::abs(lResonance.Rapidity()) > static_cast<float>(0.5))
          continue;
      }

      if (cfgCutsOnMother) {
        if (lResonance.Pt() >= cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (lResonance.M() >= cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      }

      if (cFilldeltaEtaPhiPlots) {
        if (deltaEta >= cMaxDeltaEtaCut)
          continue;
        if (deltaPhi >= cMaxDeltaPhiCut)
          continue;

        if constexpr (!IsMix) {
          histos.fill(HIST("QAafter/EtaPrafter"), trk1.eta());
          histos.fill(HIST("QAafter/PhiPrafter"), trk1.phi());
          histos.fill(HIST("QAafter/EtaKaafter"), trk2.eta());
          histos.fill(HIST("QAafter/PhiKaafter"), trk2.phi());
          histos.fill(HIST("QAafter/deltaEtaafter"), deltaEta);
          histos.fill(HIST("QAafter/deltaPhiafter"), deltaPhi);
        }
      }

      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (IsData) {
          if (isCalcRotBkg) {
            for (int i = 0; i < cNofRotations; i++) {
              float theta2 = rn->Uniform(constants::math::PI - constants::math::PI / rotationalcut, constants::math::PI + constants::math::PI / rotationalcut);
              ldaughterRot = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi() + theta2, massKa); // for rotated background
              lresonanceRot = lDecayDaughter1 + ldaughterRot;
              histos.fill(HIST("Result/Data/h3lambda1520InvMassRotation"), multiplicity, lresonanceRot.Pt(), lresonanceRot.M());
            }
          }

          if (trk1.sign() < 0) {
            if (cFillinvmass1DPlots) {
              histos.fill(HIST("Result/Data/lambda1520invmass"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmass"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (trk1.sign() > 0) {
            if (cFillinvmass1DPlots) {
              histos.fill(HIST("Result/Data/antilambda1520invmass"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3antilambda1520invmass"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        } else if (IsMix) {
          if (cFillinvmass1DPlots) {
            histos.fill(HIST("Result/Data/lambda1520invmassME"), lResonance.M());
          }
          histos.fill(HIST("Result/Data/h3lambda1520invmassME"), multiplicity, lResonance.Pt(), lResonance.M());
          if (cFilladditionalMEPlots) {
            if (trk1.sign() < 0) {
              histos.fill(HIST("Result/Data/h3lambda1520invmassME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
            } else if (trk1.sign() > 0) {
              histos.fill(HIST("Result/Data/h3lambda1520invmassME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            }
          }
        }

        // MC
        if constexpr (IsMC) {

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
          while (motherstrk1.size() > 2) {
            motherstrk1.pop_back();
            mothersPDGtrk1.pop_back();
          }

          const auto& mctrk2 = trk2.mcParticle();
          if (mctrk2.has_mothers()) {
            motherstrk2 = getMothersIndeces(mctrk2);
            mothersPDGtrk2 = getMothersPDGCodes(mctrk2);
          }
          while (motherstrk2.size() > 2) {
            motherstrk2.pop_back();
            mothersPDGtrk2.pop_back();
          }

          if (std::abs(mctrk1.pdgCode()) != 2212 || std::abs(mctrk2.pdgCode()) != 321)
            continue;

          if (motherstrk1[0] != motherstrk2[0]) // Same mother
            continue;

          if (std::abs(mothersPDGtrk1[0]) != 102134)
            continue;

          // LOGF(info, "mother trk1 id: %d, mother trk1: %d, trk1 id: %d, trk1 pdgcode: %d, mother trk2 id: %d, mother trk2: %d, trk2 id: %d, trk2 pdgcode: %d", motherstrk1[0], mothersPDGtrk1[0], trk1.globalIndex(), mctrk1.pdgCode(), motherstrk2[0], mothersPDGtrk2[0], trk2.globalIndex(), mctrk2.pdgCode());

          if (cUseEtacutMC && std::abs(lResonance.Eta()) > cEtacutMC) // eta cut
            continue;

          if (cUseRapcutMC && std::abs(lResonance.Rapidity()) > static_cast<float>(0.5)) // rapidity cut
            continue;

          histos.fill(HIST("QA/MC/h2RecoEtaPt_after"), lResonance.Eta(), lResonance.Pt());
          histos.fill(HIST("QA/MC/h2RecoPhiRapidity_after"), lResonance.Phi(), lResonance.Rapidity());

          // Track selection check.
          histos.fill(HIST("QA/MC/trkDCAxy_pr"), trk1ptPr, trk1.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAxy_ka"), trk2ptKa, trk2.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAz_pr"), trk1ptPr, trk1.dcaZ());
          histos.fill(HIST("QA/MC/trkDCAz_ka"), trk2ptKa, trk2.dcaZ());

          histos.fill(HIST("QA/MC/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          }
          histos.fill(HIST("QA/MC/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          }

          // MC histograms
          if (mothersPDGtrk1[0] > 0) {
            histos.fill(HIST("Result/MC/h3lambda1520Recoinvmass"), multiplicity, lResonance.Pt(), lResonance.M());
          } else {
            histos.fill(HIST("Result/MC/h3antilambda1520Recoinvmass"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
      } else {
        if constexpr (IsData) {
          // Like sign pair ++
          if (trk1.sign() > 0) {
            if (cFillinvmass1DPlots) {
              histos.fill(HIST("Result/Data/lambda1520invmassLSPP"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSPP"), multiplicity, lResonance.Pt(), lResonance.M());
          } else { // Like sign pair --
            if (cFillinvmass1DPlots) {
              histos.fill(HIST("Result/Data/lambda1520invmassLSMM"), lResonance.M());
            }
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSMM"), multiplicity, lResonance.Pt(), lResonance.M());
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

    colCuts.fillQA(collision);

    if (cFilladditionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);
    fillHistograms<true, false, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Lstaranalysis, processData, "Process Event for data without partition", false);

  void processMC(MCEventCandidates::iterator const& collision,
                 aod::McCollisions const&,
                 MCTrackCandidates const& tracks, aod::McParticles const&)
  {
    colCuts.fillQA(collision);

    if (std::abs(collision.posZ()) > cZvertCutMC) // Z-vertex cut
      return;

    fillHistograms<false, true, false>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Lstaranalysis, processMC, "Process Event for MC Light without partition", false);

  Partition<aod::McParticles> selectedMCParticles = (nabs(aod::mcparticle::pdgCode) == 102134); // Lambda(1520)

  void processMCTrue(MCEventCandidates::iterator const& collision, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    bool isInAfterAllCuts = colCuts.isSelected(collision);
    bool inVtx10 = (std::abs(collision.mcCollision().posZ()) > 10.) ? false : true;
    bool isTriggerTVX = collision.selection_bit(aod::evsel::kIsTriggerTVX);
    bool isSel8 = collision.sel8();
    bool isTrueINELgt0 = isTrueINEL0(collision, mcParticles);
    centrality = centEst(collision);

    auto multiplicity = collision.centFT0M();

    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());

    // Not related to the real collisions
    for (const auto& part : mcParts) { // loop over all MC particles

      if (std::abs(part.pdgCode()) != kLambda1520PDG) // Lambda1520(0)
        continue;

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

      if (cUseRapcutMC && std::abs(part.y()) > static_cast<float>(0.5)) // rapidity cut
        continue;

      histos.fill(HIST("QA/MC/h2GenEtaPt_afterRapcut"), part.eta(), part.pt());
      histos.fill(HIST("QA/MC/h2GenPhiRapidity_afterRapcut"), part.phi(), part.y());

      if (cUseEtacutMC && std::abs(part.eta()) > cEtacutMC) // eta cut
        continue;

      histos.fill(HIST("QA/MC/h2GenEtaPt_afterEtaRapCut"), part.eta(), part.pt());
      histos.fill(HIST("QA/MC/h2GenPhiRapidity_afterEtaRapCut"), part.phi(), part.y());

      // without any event selection
      if (part.pdgCode() > 0)
        histos.fill(HIST("Result/MC/Genlambda1520pt"), 0, part.pt(), multiplicity);
      else
        histos.fill(HIST("Result/MC/Genantilambda1520pt"), 0, part.pt(), multiplicity);

      if (inVtx10) // INEL10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 1, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 1, part.pt(), multiplicity);
      }
      if (inVtx10 && isSel8) // INEL>10, vtx10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 2, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 2, part.pt(), multiplicity);
      }
      if (inVtx10 && isTriggerTVX) // vtx10, TriggerTVX
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 3, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 3, part.pt(), multiplicity);
      }
      if (isInAfterAllCuts) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("Result/MC/Genlambda1520pt"), 4, part.pt(), multiplicity);
        else
          histos.fill(HIST("Result/MC/Genantilambda1520pt"), 4, part.pt(), multiplicity);
      }
    }

    // QA for Trigger efficiency
    histos.fill(HIST("Event/hMCEventIndices"), centrality, kINEL);
    if (inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kINEL10);
    if (isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kINELg0);
    if (inVtx10 && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kINELg010);

    // TVX MB trigger
    if (isTriggerTVX)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kTrig);
    if (isTriggerTVX && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kTrig10);
    if (isTriggerTVX && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kTrigINELg0);
    if (isTriggerTVX && isTrueINELgt0 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kTrigINELg010);

    // Sel8 event selection
    if (isSel8)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kSel8);
    if (isSel8 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kSel810);
    if (isSel8 && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kSel8INELg0);
    if (isSel8 && isTrueINELgt0 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kSel8INELg010);

    // CollisionCuts selection
    if (isInAfterAllCuts)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kAllCuts);
    if (isInAfterAllCuts && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kAllCuts10);
    if (isInAfterAllCuts && isTrueINELgt0)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kAllCutsINELg0);
    if (isInAfterAllCuts && isTrueINELgt0 && inVtx10)
      histos.fill(HIST("Event/hMCEventIndices"), centrality, kAllCutsINELg010);
  }
  PROCESS_SWITCH(Lstaranalysis, processMCTrue, "Process Event for MC only", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<collision::PosZ, cent::CentFT0M>;
  BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};

  void processME(EventCandidates const& collision,
                 TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collision, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      // LOGF(info, "Mixed event collisions: (%d, %d)", collision1.globalIndex(), collision2.globalIndex());

      // for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
      //  LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.index(), t2.index(), collision1.index(), collision2.index());
      //  }

      if (cFilladditionalQAeventPlots) {
        // Fill histograms for the characteristics of the *mixed* events (collision1 and collision2)
        // This will show the distribution of events that are actually being mixed.
        histos.fill(HIST("QAevent/hMixPool_VtxZ"), collision1.posZ());
        histos.fill(HIST("QAevent/hMixPool_Multiplicity"), collision1.centFT0M()); // Assuming getCentrality() gives multiplicity
        histos.fill(HIST("QAevent/hMixPool_VtxZ_vs_Multiplicity"), collision1.posZ(), collision1.centFT0M());

        // You might also want to fill for collision2 if you want to see both partners' distributions
        // histos.fill(HIST("QAevent/hMixPool_VtxZ"), collision2.posZ());
        // histos.fill(HIST("QAevent/hMixPool_Multiplicity"), collision2.getCentrality());
        // histos.fill(HIST("QAevent/hMixPool_VtxZ_vs_Multiplicity"), collision2.posZ(), collision2.getCentrality());

        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.f);
      }
      fillHistograms<false, false, true>(collision1, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(Lstaranalysis, processME, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Lstaranalysis>(cfgc)};
}
