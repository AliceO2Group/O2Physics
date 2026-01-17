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
/// \file jFlucEfficiencyTask.cxx
/// \brief Task to calculate the efficiency of the cf-derived tracks/particles
/// \author DongJo Kim, Jasper Parkkila, Bong-Hwi Lim (djkim@cern.ch, jparkkil@cern.ch, bong-hwi.lim@cern.ch)
/// \since March 2024

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGLF/Utils/collisionCuts.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::rctsel;
using namespace o2::constants::physics;

struct JFlucEfficiencyTask {
  Service<o2::framework::O2DatabasePDG> pdg;
  // Add the pT binning array as a static member
  static constexpr std::array<double, 94> PttJacek = {
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
    1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
    4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
    26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0,
    70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0,
    170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0,
    270.0, 280.0, 290.0, 300.0};

  static constexpr double kChargeThreshold = 3.0; // PDG charge units: 3 = |e|

  // Update the axisPt configuration with proper vector initialization
  ConfigurableAxis axisPt{"axisPt", std::vector<double>(PttJacek.begin(), PttJacek.end()), "pT axis"};

  // Event cuts
  Configurable<int> cfgAcceptSplitCollisions{"cfgAcceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  o2::analysis::CollisonCuts colCuts;
  struct : ConfigurableGroup {
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Min centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 100.0f, "Max centrality"};
    Configurable<int> cfgEvtOccupancyInTimeRangeMax{"cfgEvtOccupancyInTimeRangeMax", -1, "Evt sel: maximum track occupancy"};
    Configurable<int> cfgEvtOccupancyInTimeRangeMin{"cfgEvtOccupancyInTimeRangeMin", -1, "Evt sel: minimum track occupancy"};
    Configurable<bool> cfgEvtTriggerCheck{"cfgEvtTriggerCheck", false, "Evt sel: check for trigger"};
    Configurable<bool> cfgEvtOfflineCheck{"cfgEvtOfflineCheck", true, "Evt sel: check for offline selection"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", false, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtTFBorderCut{"cfgEvtTFBorderCut", false, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtUseITSTPCvertex{"cfgEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
    Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", true, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", true, "Evt sel: apply pileup rejection"};
    Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtCollInTimeRangeStandard{"cfgEvtCollInTimeRangeStandard", true, "Evt sel: apply NoCollInTimeRangeStandard"};
    Configurable<bool> cfgEvtRun2AliEventCuts{"cfgEvtRun2AliEventCuts", true, "Evt sel: apply Run2 Ali event cuts"};
    Configurable<bool> cfgEvtRun2INELgtZERO{"cfgEvtRun2INELgtZERO", false, "Evt sel: apply Run2 INEL>0 event cuts"};
    Configurable<bool> cfgEvtUseRCTFlagChecker{"cfgEvtUseRCTFlagChecker", false, "Evt sel: use RCT flag checker"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } EventCuts;
  RCTFlagsChecker rctChecker;

  // Track selections
  struct : ConfigurableGroup {
    Configurable<float> cfgMinPt{"cfgMinPt", 0.6, "Track minium pt cut"};
    Configurable<float> cfgMaxPt{"cfgMaxPt", 300.0f, "Maximum transverse momentum"};
    Configurable<float> cfgEtaMin{"cfgEtaMin", -1.0f, "Minimum pseudorapidity"};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 1.0f, "Maximum pseudorapidity"};
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", true, "Global track selection"};                        // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};           // PV Contriuibutor
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
    Configurable<float> cfgMaxbDCArToPVcut{"cfgMaxbDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
    Configurable<float> cfgMaxbDCAzToPVcut{"cfgMaxbDCAzToPVcut", 1.0, "Track DCAz cut to PV Maximum"};
    // PID
    Configurable<float> cfgPIDnSigmaCut{"cfgPIDnSigmaCut", 3.0, "PID nSigma cut"};
  } TrackCuts;

  Configurable<bool> applyMCStudy{"applyMCStudy", false, "Apply MC study"};

  // Configurable for track selection
  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

  // Configurable axes
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "multiplicity / centrality axis"};

  // Filter declarations
  Filter cfCollisionFilter = (nabs(aod::collision::posZ) < EventCuts.cfgEvtZvtx);
  Filter cfTrackFilter = (aod::cftrack::pt >= TrackCuts.cfgMinPt) &&
                         (aod::cftrack::pt <= TrackCuts.cfgMaxPt) &&
                         (aod::cftrack::eta >= TrackCuts.cfgEtaMin) &&
                         (aod::cftrack::eta <= TrackCuts.cfgEtaMax);
  // Filter collisionFilter = (nabs(aod::collision::posZ) < EventCuts.cfgEvtZvtx);
  Filter trackFilter = (aod::track::pt >= TrackCuts.cfgMinPt) &&
                       (aod::track::pt <= TrackCuts.cfgMaxPt) &&
                       (aod::track::eta >= TrackCuts.cfgEtaMin) &&
                       (aod::track::eta <= TrackCuts.cfgEtaMax);

  Configurable<int> cfgCentBinsForMC{"cfgCentBinsForMC", 1, "Centrality bins for MC, 0: off, 1: on"};
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As>;
  using CollisionRun2Candidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksCov>;
  using TrackCandidatesPID = soa::Join<TrackCandidates, aod::pidTPCPi>;
  using MCCollisionCandidates = soa::Join<CollisionCandidates, aod::McCollisionLabel>;
  using MCRun2CollisionCandidates = soa::Join<CollisionRun2Candidates, aod::McCollisionLabel>;
  using MCTrackCandidates = soa::Join<TrackCandidates, aod::McTrackLabel>;
  using MCTrackCandidatesPID = soa::Join<MCTrackCandidates, aod::pidTPCPi>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  // Histogram Registry
  HistogramRegistry registry{
    "registry",
    {{"hEventCounterMC", "Event counter MC;Counter;Counts", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hEventCounterReco", "Event counter Reco;Counter;Counts", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hZVertexMC", "MC Z vertex distribution;Z vertex (cm);Centrality (%)", {HistType::kTH2F, {{200, -20, 20}, {axisMultiplicity}}}},
     {"hZVertexReco", "Reconstructed Z vertex distribution;Z vertex (cm);Centrality (%)", {HistType::kTH2F, {{200, -20, 20}, {axisMultiplicity}}}},
     {"hZVertexCorrelation", "Z vertex correlation;MC Z vertex (cm);Reco Z vertex (cm)", {HistType::kTH2F, {{200, -20, 20}, {200, -20, 20}}}}}};

  // Configurable for debugging
  Configurable<bool> debugMode{"debugMode", false, "Debug mode"};

  void init(InitContext const&)
  {
    if (debugMode) {
      LOGF(info, "Initializing JFlucEfficiencyTask");
    }
    if (!doprocessMCRun2 && !doprocessDataRun2) {
      colCuts.setCuts(EventCuts.cfgEvtZvtx, EventCuts.cfgEvtTriggerCheck, EventCuts.cfgEvtOfflineCheck, /*checkRun3*/ true, /*triggerTVXsel*/ false, EventCuts.cfgEvtOccupancyInTimeRangeMax, EventCuts.cfgEvtOccupancyInTimeRangeMin);
    } else {
      colCuts.setCuts(EventCuts.cfgEvtZvtx, EventCuts.cfgEvtTriggerCheck, EventCuts.cfgEvtOfflineCheck, false);
    }
    colCuts.init(&registry);
    colCuts.setTriggerTVX(EventCuts.cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(EventCuts.cfgEvtTFBorderCut);
    colCuts.setApplyITSTPCvertex(EventCuts.cfgEvtUseITSTPCvertex);
    colCuts.setApplyZvertexTimedifference(EventCuts.cfgEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(EventCuts.cfgEvtPileupRejection);
    colCuts.setApplyNoITSROBorderCut(EventCuts.cfgEvtNoITSROBorderCut);
    colCuts.setApplyCollInTimeRangeStandard(EventCuts.cfgEvtCollInTimeRangeStandard);
    colCuts.setApplyRun2AliEventCuts(EventCuts.cfgEvtRun2AliEventCuts);
    colCuts.setApplyRun2INELgtZERO(EventCuts.cfgEvtRun2INELgtZERO);
    colCuts.printCuts();

    rctChecker.init(EventCuts.cfgEvtRCTFlagCheckerLabel, EventCuts.cfgEvtRCTFlagCheckerZDCCheck, EventCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    // Helper function to add histograms with consistent naming
    auto addHistograms = [this](const std::string& prefix, bool isMC = false) {
      if (isMC) {
        // Generated (MC) histograms - pT has all variations
        registry.add(Form("hPtGen%s", prefix.c_str()),
                     Form("Generated p_{T} %s;p_{T} (GeV/c);Centrality (%%);Counts", prefix.c_str()),
                     {HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      }

      // Reconstructed histograms - pT has all variations
      registry.add(Form("hPtRec%s", prefix.c_str()),
                   Form("Reconstructed p_{T} %s;p_{T} (GeV/c);Centrality (%%);Counts", prefix.c_str()),
                   {HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
    };

    // Add MC histograms if MC processing is enabled
    if (doprocessDerivedMC || doprocessMC || doprocessMCRun2 || doprocessMCPID) {
      addHistograms("", true);        // hPtGen, hPtRec
      addHistograms("Pos", true);     // hPtGenPos, hPtRecPos
      addHistograms("Neg", true);     // hPtGenNeg, hPtRecNeg
      addHistograms("_Pos", true);    // hPtGen_Pos, hPtRec_Pos
      addHistograms("_Neg", true);    // hPtGen_Neg, hPtRec_Neg
      addHistograms("Pos_Pos", true); // hPtGenPos_Pos, hPtRecPos_Pos
      addHistograms("Pos_Neg", true); // hPtGenPos_Neg, hPtRecPos_Neg
      addHistograms("Neg_Pos", true); // hPtGenNeg_Pos, hPtRecNeg_Pos
      addHistograms("Neg_Neg", true); // hPtGenNeg_Neg, hPtRecNeg_Neg
    } else {
      // Add reconstructed histograms
      addHistograms("");        // hPtRec
      addHistograms("_Pos");    // hPtRec_Pos
      addHistograms("_Neg");    // hPtRec_Neg
      addHistograms("Pos");     // hPtRecPos
      addHistograms("Neg");     // hPtRecNeg
      addHistograms("Pos_Pos"); // hPtRecPos_Pos
      addHistograms("Pos_Neg"); // hPtRecPos_Neg
      addHistograms("Neg_Pos"); // hPtRecNeg_Pos
      addHistograms("Neg_Neg"); // hPtRecNeg_Neg
    }
    // Add basic eta histograms separately
    if (doprocessDerivedMC || doprocessMC || doprocessMCRun2 || doprocessMCPID) {
      registry.add("hEtaGen", "Generated #eta (all);#eta;Centrality (%);Counts",
                   {HistType::kTH2F, {AxisSpec(100, -1, 1), AxisSpec(axisMultiplicity)}});
    }
    registry.add("hEtaRec", "Reconstructed #eta (all);#eta;Centrality (%);Counts",
                 {HistType::kTH2F, {AxisSpec(100, -1, 1), AxisSpec(axisMultiplicity)}});

    // Add pT uncertainty histograms
    registry.add("hPtUncertainty", "Track p_{T} uncertainty;p_{T} (GeV/c);Centrality (%);p_{T} uncertainty (GeV/c)",
                 {HistType::kTH3F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity), AxisSpec(100, 0, 10)}});
    registry.add("hPtUncertaintyPos", "Track p_{T} uncertainty (positive);p_{T} (GeV/c);Centrality (%);p_{T} uncertainty (GeV/c)",
                 {HistType::kTH3F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity), AxisSpec(100, 0, 10)}});
    registry.add("hPtUncertaintyNeg", "Track p_{T} uncertainty (negative);p_{T} (GeV/c);Centrality (%);p_{T} uncertainty (GeV/c)",
                 {HistType::kTH3F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity), AxisSpec(100, 0, 10)}});

    // Add MC study histograms if enabled
    if (applyMCStudy && (doprocessDerivedMC || doprocessMC || doprocessMCRun2 || doprocessMCPID)) {
      registry.add("hChargeSignMismatch", "Charge-Sign mismatch cases", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMismatchPos", "MC charge + but track sign -", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMismatchNeg", "MC charge - but track sign +", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMatch", "Charge-Sign match cases", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMatchPos", "MC charge + and track sign +", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMatchNeg", "MC charge - and track sign -", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
    }

    // Add efficiency histograms if enabled
    if (doprocessEfficiency) {
      registry.add("hPtGenData", "Generated p_{T} from data events (all);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {axisPt, axisMultiplicity}});
      registry.add("hEtaGenData", "Generated #eta from data events (all);#eta;Centrality (%);Counts",
                   {HistType::kTH2F, {AxisSpec(100, -1, 1), axisMultiplicity}});
      registry.add("hPtGenDataPos", "Generated p_{T} from data events (positive);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {axisPt, axisMultiplicity}});
      registry.add("hPtGenDataNeg", "Generated p_{T} from data events (negative);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {axisPt, axisMultiplicity}});
      registry.add("hPtRecData", "Reconstructed p_{T} (all);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hEtaRecData", "Reconstructed #eta (all);#eta;Centrality (%);Counts",
                   {HistType::kTH2F, {AxisSpec(100, -1, 1), AxisSpec(axisMultiplicity)}});
      registry.add("hPtRecDataPos", "Reconstructed p_{T} (positive);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hPtRecDataNeg", "Reconstructed p_{T} (negative);p_{T} (GeV/c);Centrality (%);Counts",
                   {HistType::kTH2F, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
    }

    // Histogram labels
    auto h1 = registry.get<TH1>(HIST("hEventCounterMC"));
    auto h2 = registry.get<TH1>(HIST("hEventCounterReco"));

    if (h1 && h2) {
      h1->GetXaxis()->SetBinLabel(1, "All MC Events");
      h1->GetXaxis()->SetBinLabel(2, "Selected MC Events");
      h1->GetXaxis()->SetBinLabel(3, "Analyzed MC Events");

      h2->GetXaxis()->SetBinLabel(1, "All Reco Events");
      h2->GetXaxis()->SetBinLabel(2, "Selected Reco Events");
      h2->GetXaxis()->SetBinLabel(3, "Analyzed Reco Events");
    } else {
      LOGF(error, "Failed to get histograms from registry");
    }
  }

  template <typename ParticleType>
  double getCharge(ParticleType const& particle)
  {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle) {
      return 10.f;
    }
    return pdgParticle->Charge();
  }
  bool isChargedParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= kChargeThreshold;
  }
  // Track selection
  template <typename TrackType>
  bool trackCut(TrackType const& track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < TrackCuts.cfgMinPt)
      return false;
    if (std::abs(track.pt()) > TrackCuts.cfgMaxPt)
      return false;
    if (track.eta() < TrackCuts.cfgEtaMin)
      return false;
    if (track.eta() > TrackCuts.cfgEtaMax)
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

  void processDerivedMC(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Filtered<aod::CFMcParticles> const& mcParticles)
  {
    float centrality = mcCollision.multiplicity(); // multiplicity: number of primary particles TODO: apply percentiles
    registry.fill(HIST("hEventCounterMC"), 0);
    registry.fill(HIST("hZVertexMC"), mcCollision.posZ(), centrality);

    for (const auto& particle : mcParticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      registry.fill(HIST("hPtGen"), particle.pt(), centrality);
      registry.fill(HIST("hEtaGen"), particle.eta(), centrality);

      if (particle.sign() > 0) { // Positive particles
        registry.fill(HIST("hPtGenPos"), particle.pt(), centrality);
      } else if (particle.sign() < 0) { // Negative particles
        registry.fill(HIST("hPtGenNeg"), particle.pt(), centrality);
      }
    }
  }

  void processDerivedData(soa::Filtered<aod::CFCollisions>::iterator const& cfCollision, soa::Filtered<aod::CFTracks> const& cfTracks)
  {
    float centrality = cfCollision.multiplicity();

    if (centrality < EventCuts.cfgCentMin || centrality > EventCuts.cfgCentMax) {
      return;
    }
    registry.fill(HIST("hZVertexReco"), cfCollision.posZ(), centrality);

    for (const auto& track : cfTracks) {
      registry.fill(HIST("hPtRec"), track.pt(), centrality);
      registry.fill(HIST("hEtaRec"), track.eta(), centrality);

      if (track.sign() > 0) { // Positive tracks
        registry.fill(HIST("hPtRecPos"), track.pt(), centrality);
      } else if (track.sign() < 0) { // Negative tracks
        registry.fill(HIST("hPtRecNeg"), track.pt(), centrality);
      }
    }
  }

  Preslice<TrackCandidates> perCollision = aod::track::collisionId;
  // Common histogram filling function for tracks
  template <typename TrackType>
  void fillTrackHistograms(const TrackType& track, float centrality)
  {
    // Basic pT and eta histograms
    registry.fill(HIST("hPtRec"), track.pt(), centrality);
    registry.fill(HIST("hEtaRec"), track.eta(), centrality);

    // pT histograms by eta direction
    if (track.eta() > 0) {
      registry.fill(HIST("hPtRec_Pos"), track.pt(), centrality);
    } else if (track.eta() < 0) {
      registry.fill(HIST("hPtRec_Neg"), track.pt(), centrality);
    }

    // Charge sign specific histograms
    if (track.sign() > 0) { // Positive tracks
      registry.fill(HIST("hPtRecPos"), track.pt(), centrality);
      if (track.eta() > 0) {
        registry.fill(HIST("hPtRecPos_Pos"), track.pt(), centrality);
      } else if (track.eta() < 0) {
        registry.fill(HIST("hPtRecPos_Neg"), track.pt(), centrality);
      }
    } else if (track.sign() < 0) { // Negative tracks
      registry.fill(HIST("hPtRecNeg"), track.pt(), centrality);
      if (track.eta() > 0) {
        registry.fill(HIST("hPtRecNeg_Pos"), track.pt(), centrality);
      } else if (track.eta() < 0) {
        registry.fill(HIST("hPtRecNeg_Neg"), track.pt(), centrality);
      }
    }

    // pT uncertainty histograms
    auto ptUncertaintySigma = track.c1Pt21Pt2() * track.pt() * track.pt(); // Variance of pT
    auto ptUncertainty = std::sqrt(ptUncertaintySigma);                    // Standard deviation of pT
    registry.fill(HIST("hPtUncertainty"), track.pt(), centrality, ptUncertainty);

    if (track.sign() > 0) {
      registry.fill(HIST("hPtUncertaintyPos"), track.pt(), centrality, ptUncertainty);
    } else if (track.sign() < 0) {
      registry.fill(HIST("hPtUncertaintyNeg"), track.pt(), centrality, ptUncertainty);
    }
  }

  // Common histogram filling function for MC particles
  template <typename ParticleType>
  void fillMCParticleHistograms(const ParticleType& particle, float centrality)
  {
    registry.fill(HIST("hPtGen"), particle.pt(), centrality);
    registry.fill(HIST("hEtaGen"), particle.eta(), centrality);

    // pT histograms by eta direction
    if (particle.eta() > 0) {
      registry.fill(HIST("hPtGen_Pos"), particle.pt(), centrality);
    } else if (particle.eta() < 0) {
      registry.fill(HIST("hPtGen_Neg"), particle.pt(), centrality);
    }

    // Charge sign specific histograms
    auto charge = getCharge(particle);
    if (charge > 0) { // Positive particles
      registry.fill(HIST("hPtGenPos"), particle.pt(), centrality);
      if (particle.eta() > 0) {
        registry.fill(HIST("hPtGenPos_Pos"), particle.pt(), centrality);
      } else if (particle.eta() < 0) {
        registry.fill(HIST("hPtGenPos_Neg"), particle.pt(), centrality);
      }
    } else if (charge < 0) { // Negative particles
      registry.fill(HIST("hPtGenNeg"), particle.pt(), centrality);
      if (particle.eta() > 0) {
        registry.fill(HIST("hPtGenNeg_Pos"), particle.pt(), centrality);
      } else if (particle.eta() < 0) {
        registry.fill(HIST("hPtGenNeg_Neg"), particle.pt(), centrality);
      }
    }
  }

  // Common event selection function
  template <typename CollisionType>
  bool selectEvent(const CollisionType& collision, float& centrality)
  {
    if (!colCuts.isSelected(collision)) {
      return false;
    }

    if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }

    colCuts.fillQA(collision);
    centrality = collision.centFT0C();

    if (centrality < EventCuts.cfgCentMin || centrality > EventCuts.cfgCentMax) {
      return false;
    }

    return true;
  }

  // Common event selection function for Run2 data
  template <typename CollisionType>
  bool selectEventRun2(const CollisionType& collision, float& centrality)
  {
    if (!colCuts.isSelected(collision)) {
      return false;
    }

    colCuts.fillQARun2(collision);
    centrality = collision.centRun2V0M();

    if (centrality < EventCuts.cfgCentMin || centrality > EventCuts.cfgCentMax) {
      return false;
    }

    return true;
  }

  void processMC(aod::McCollisions::iterator const& mcCollision,
                 soa::SmallGroups<MCCollisionCandidates> const& collisions,
                 soa::Filtered<MCTrackCandidates> const& mcTracks,
                 aod::McParticles const& mcParticles)
  {
    registry.fill(HIST("hEventCounterMC"), 0);
    if (!(std::abs(mcCollision.posZ()) < EventCuts.cfgEvtZvtx)) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    if (cfgAcceptSplitCollisions == 0 && collisions.size() > 1) {
      return;
    }

    float centrality = -999;
    for (const auto& collision : collisions) { // Anyway only 1 collision per mcCollision will be selected
      if (!colCuts.isSelected(collision))      // Default event selection
        return;
      if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(collision)) {
        return;
      }
      colCuts.fillQA(collision);
      centrality = collision.centFT0C();
    }

    registry.fill(HIST("hEventCounterMC"), 1);
    registry.fill(HIST("hZVertexMC"), mcCollision.posZ(), centrality);

    if (centrality < EventCuts.cfgCentMin || centrality > EventCuts.cfgCentMax) {
      return;
    }

    // Fill MC particle histograms
    for (const auto& particle : mcParticles) {
      if ((!particle.isPhysicalPrimary()) || !isChargedParticle(particle.pdgCode())) {
        continue;
      }
      // pT and eta selections
      if (particle.pt() < TrackCuts.cfgMinPt || particle.pt() > TrackCuts.cfgMaxPt ||
          particle.eta() < TrackCuts.cfgEtaMin || particle.eta() > TrackCuts.cfgEtaMax) {
        continue;
      }
      fillMCParticleHistograms(particle, centrality);
    }

    // Process reconstructed tracks
    for (const auto& collision : collisions) {
      registry.fill(HIST("hZVertexReco"), collision.posZ(), centrality);
      registry.fill(HIST("hZVertexCorrelation"), mcCollision.posZ(), collision.posZ());

      auto tracks = mcTracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        if (!trackCut(track)) {
          continue;
        }
        auto mcPart = track.mcParticle();
        if (!mcPart.isPhysicalPrimary() || !isChargedParticle(mcPart.pdgCode())) {
          continue;
        }

        if (applyMCStudy) {
          // Check charge-sign consistency
          auto mcCharge = getCharge(mcPart);
          auto trackSign = track.sign();

          if (mcCharge > 0 && trackSign > 0) {
            // MC charge + and track sign +
            registry.fill(HIST("hChargeSignMatchPos"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMatch"), track.pt(), centrality);
          } else if (mcCharge < 0 && trackSign < 0) {
            // MC charge - and track sign -
            registry.fill(HIST("hChargeSignMatchNeg"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMatch"), track.pt(), centrality);
          } else if (mcCharge > 0 && trackSign < 0) {
            // MC charge + but track sign -
            registry.fill(HIST("hChargeSignMismatchPos"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMismatch"), track.pt(), centrality);
          } else if (mcCharge < 0 && trackSign > 0) {
            // MC charge - but track sign +
            registry.fill(HIST("hChargeSignMismatchNeg"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMismatch"), track.pt(), centrality);
          }
        }

        fillTrackHistograms(track, centrality);
      }
    }
  }

  Preslice<TrackCandidatesPID> perCollisionPID = aod::track::collisionId;
  void processMCPID(aod::McCollisions::iterator const& mcCollision,
                    soa::SmallGroups<MCCollisionCandidates> const& collisions,
                    soa::Filtered<MCTrackCandidatesPID> const& mcTracks,
                    aod::McParticles const& mcParticles)
  {
    registry.fill(HIST("hEventCounterMC"), 0);
    if (!(std::abs(mcCollision.posZ()) < EventCuts.cfgEvtZvtx)) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    if (cfgAcceptSplitCollisions == 0 && collisions.size() > 1) {
      return;
    }

    float centrality = -999;
    for (const auto& collision : collisions) { // Anyway only 1 collision per mcCollision will be selected
      if (!colCuts.isSelected(collision))      // Default event selection
        return;
      if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(collision)) {
        return;
      }
      colCuts.fillQA(collision);
      centrality = collision.centFT0C();
    }

    registry.fill(HIST("hEventCounterMC"), 1);
    registry.fill(HIST("hZVertexMC"), mcCollision.posZ(), centrality);

    if (centrality < EventCuts.cfgCentMin || centrality > EventCuts.cfgCentMax) {
      return;
    }

    // Fill MC particle histograms
    for (const auto& particle : mcParticles) {
      if (!isChargedParticle(particle.pdgCode())) {
        continue;
      }
      // pT and eta selections
      if (particle.pt() < TrackCuts.cfgMinPt || particle.pt() > TrackCuts.cfgMaxPt ||
          particle.eta() < TrackCuts.cfgEtaMin || particle.eta() > TrackCuts.cfgEtaMax) {
        continue;
      }
      // Check if the particle is a pion
      if (std::abs(particle.pdgCode()) != std::abs(kPiPlus)) {
        continue;
      }
      fillMCParticleHistograms(particle, centrality);
    }

    // Process reconstructed tracks
    for (const auto& collision : collisions) {
      registry.fill(HIST("hZVertexReco"), collision.posZ(), centrality);
      registry.fill(HIST("hZVertexCorrelation"), mcCollision.posZ(), collision.posZ());

      auto tracks = mcTracks.sliceBy(perCollisionPID, collision.globalIndex());
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        if (std::abs(track.tpcNSigmaPi()) > TrackCuts.cfgPIDnSigmaCut) {
          continue;
        }
        if (!trackCut(track)) {
          continue;
        }
        auto mcPart = track.mcParticle();
        if (!isChargedParticle(mcPart.pdgCode())) {
          continue;
        }
        if (std::abs(mcPart.pdgCode()) != kPiPlus) {
          continue;
        }

        if (applyMCStudy) {
          // Check charge-sign consistency
          auto mcCharge = getCharge(mcPart);
          auto trackSign = track.sign();

          if (mcCharge > 0 && trackSign > 0) {
            // MC charge + and track sign +
            registry.fill(HIST("hChargeSignMatchPos"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMatch"), track.pt(), centrality);
          } else if (mcCharge < 0 && trackSign < 0) {
            // MC charge - and track sign -
            registry.fill(HIST("hChargeSignMatchNeg"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMatch"), track.pt(), centrality);
          } else if (mcCharge > 0 && trackSign < 0) {
            // MC charge + but track sign -
            registry.fill(HIST("hChargeSignMismatchPos"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMismatch"), track.pt(), centrality);
          } else if (mcCharge < 0 && trackSign > 0) {
            // MC charge - but track sign +
            registry.fill(HIST("hChargeSignMismatchNeg"), track.pt(), centrality);
            registry.fill(HIST("hChargeSignMismatch"), track.pt(), centrality);
          }
        }

        fillTrackHistograms(track, centrality);
      }
    }
  }

  void processMCRun2(aod::McCollisions::iterator const& mcCollision,
                     soa::SmallGroups<MCRun2CollisionCandidates> const& collisions,
                     soa::Filtered<MCTrackCandidates> const& mcTracks,
                     aod::McParticles const& mcParticles,
                     BCsWithRun2Info const&)
  {
    registry.fill(HIST("hEventCounterMC"), 0);
    if (!(std::abs(mcCollision.posZ()) < EventCuts.cfgEvtZvtx)) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    if (cfgAcceptSplitCollisions == 0 && collisions.size() > 1) {
      return;
    }

    float centrality = -999;
    for (const auto& collision : collisions) { // Anyway only 1 collision per mcCollision will be selected
      if (!colCuts.isSelected(collision))      // Default event selection
        return;
      colCuts.fillQARun2(collision);
      centrality = collision.centRun2V0M();
    }

    registry.fill(HIST("hEventCounterMC"), 1);
    registry.fill(HIST("hZVertexMC"), mcCollision.posZ(), centrality);

    if (centrality < EventCuts.cfgCentMin || centrality > EventCuts.cfgCentMax) {
      return;
    }

    // Fill MC particle histograms
    for (const auto& particle : mcParticles) {
      if ((!particle.isPhysicalPrimary()) || !isChargedParticle(particle.pdgCode())) {
        continue;
      }
      // pT and eta selections
      if (particle.pt() < TrackCuts.cfgMinPt || particle.pt() > TrackCuts.cfgMaxPt ||
          particle.eta() < TrackCuts.cfgEtaMin || particle.eta() > TrackCuts.cfgEtaMax) {
        continue;
      }
      fillMCParticleHistograms(particle, centrality);
    }

    // Process reconstructed tracks
    for (const auto& collision : collisions) {
      registry.fill(HIST("hZVertexReco"), collision.posZ(), centrality);
      registry.fill(HIST("hZVertexCorrelation"), mcCollision.posZ(), collision.posZ());

      auto tracks = mcTracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto mcPart = track.mcParticle();
        if (!mcPart.isPhysicalPrimary() || !isChargedParticle(mcPart.pdgCode())) {
          continue;
        }
        if (!trackCut(track)) {
          continue;
        }

        fillTrackHistograms(track, centrality);
      }
    }
  }

  void processData(CollisionCandidates::iterator const& collision, soa::Filtered<TrackCandidates> const& tracks)
  {
    float centrality;
    if (!selectEvent(collision, centrality)) {
      return;
    }

    registry.fill(HIST("hZVertexReco"), collision.posZ(), centrality);

    for (const auto& track : tracks) {
      if (!trackCut(track)) {
        continue;
      }
      fillTrackHistograms(track, centrality);
    }
  }

  void processDataPID(CollisionCandidates::iterator const& collision, soa::Filtered<TrackCandidatesPID> const& tracks)
  {
    float centrality;
    if (!selectEvent(collision, centrality)) {
      return;
    }

    registry.fill(HIST("hZVertexReco"), collision.posZ(), centrality);

    for (const auto& track : tracks) {
      if (!trackCut(track)) {
        continue;
      }
      if (std::abs(track.tpcNSigmaPi()) > TrackCuts.cfgPIDnSigmaCut) {
        continue;
      }
      fillTrackHistograms(track, centrality);
    }
  }

  void processDataRun2(CollisionRun2Candidates::iterator const& collision, soa::Filtered<TrackCandidates> const& tracks, BCsWithRun2Info const&)
  {
    float centrality;
    if (!selectEventRun2(collision, centrality)) {
      return;
    }

    registry.fill(HIST("hZVertexReco"), collision.posZ(), centrality);

    for (const auto& track : tracks) {
      if (!trackCut(track)) {
        continue;
      }
      fillTrackHistograms(track, centrality);
    }
  }

  // NOTE SmallGroups includes soa::Filtered always
  Preslice<aod::CFTracksWithLabel> perCFCollision = aod::cftrack::cfCollisionId;
  void processEfficiency(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision,
                         aod::CFMcParticles const& mcParticles,
                         soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions,
                         soa::Filtered<aod::CFTracksWithLabel> const& tracks)
  {
    try {
      // Count MC events and fill MC z-vertex with centrality
      if (debugMode) {
        LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
      }
      auto multiplicity = mcCollision.multiplicity();
      if (cfgCentBinsForMC > 0) {
        if (collisions.size() == 0) {
          return;
        }
        for (const auto& collision : collisions) {
          multiplicity = collision.multiplicity();
        }
      }
      if (debugMode) {
        LOGF(info, "MC collision multiplicity: %f", multiplicity);
      }
      registry.fill(HIST("hEventCounterMC"), 0);
      registry.fill(HIST("hZVertexMC"), mcCollision.posZ(), multiplicity);
      if (debugMode) {
        LOGF(info, "Processing MC collision %d at z = %.3f", mcCollision.globalIndex(), mcCollision.posZ());
      }

      // Fill MC particle histograms
      for (const auto& mcParticle : mcParticles) {
        if (!mcParticle.isPhysicalPrimary() || mcParticle.sign() == 0) {
          continue;
        }
        registry.fill(HIST("hPtGenData"), mcParticle.pt(), multiplicity);
        registry.fill(HIST("hEtaGenData"), mcParticle.eta(), multiplicity);
        if (mcParticle.sign() > 0) {
          registry.fill(HIST("hPtGenDataPos"), mcParticle.pt(), multiplicity);
        } else if (mcParticle.sign() < 0) {
          registry.fill(HIST("hPtGenDataNeg"), mcParticle.pt(), multiplicity);
        }
      }
      registry.fill(HIST("hEventCounterMC"), 1);

      // Check reconstructed collisions
      if (collisions.size() == 0) {
        if (debugMode) {
          LOGF(info, "No reconstructed collisions found for MC collision %d", mcCollision.globalIndex());
        }
        return;
      }

      // Process reconstructed events
      for (const auto& collision : collisions) {
        registry.fill(HIST("hEventCounterReco"), 0);
        registry.fill(HIST("hZVertexReco"), collision.posZ(), collision.multiplicity());
        registry.fill(HIST("hZVertexCorrelation"), mcCollision.posZ(), collision.posZ());

        if (debugMode) {
          LOGF(info, "Processing reconstructed collision %d at z = %.3f",
               collision.globalIndex(), collision.posZ());
        }
        registry.fill(HIST("hEventCounterReco"), 1);

        // Fill track histograms
        auto groupedTracks = tracks.sliceBy(perCFCollision, collision.globalIndex());
        if (debugMode) {
          LOGF(info, "Reconstructed collision %d has %d tracks", collision.globalIndex(), groupedTracks.size());
        }
        for (const auto& track : groupedTracks) {
          if (!track.has_cfMCParticle()) {
            if (debugMode) {
              LOGF(debug, "Track without MC particle found");
            }
            continue;
          }
          // primary particles only
          const auto& mcParticle = track.cfMCParticle();
          if (!mcParticle.isPhysicalPrimary()) {
            continue;
          }
          registry.fill(HIST("hPtRecData"), track.pt(), collision.multiplicity());
          registry.fill(HIST("hEtaRecData"), track.eta(), collision.multiplicity());
          if (track.sign() > 0) {
            registry.fill(HIST("hPtRecDataPos"), track.pt(), collision.multiplicity());
          } else if (track.sign() < 0) {
            registry.fill(HIST("hPtRecDataNeg"), track.pt(), collision.multiplicity());
          }
        }

        // Count selected and analyzed events
        registry.fill(HIST("hEventCounterReco"), 2);
      }

      registry.fill(HIST("hEventCounterMC"), 2);

    } catch (const std::exception& e) {
      LOGF(error, "Exception caught in processEfficiency: %s", e.what());
    } catch (...) {
      LOGF(error, "Unknown exception caught in processEfficiency");
    }
  }

  PROCESS_SWITCH(JFlucEfficiencyTask, processMC, "Process MC only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processMCPID, "Process MC with PID only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processMCRun2, "Process Run2 MC only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processData, "Process data only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processDataPID, "Process data with PID only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processDataRun2, "Process Run2 data only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processDerivedMC, "Process derived MC only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processDerivedData, "Process derived data only", false);
  PROCESS_SWITCH(JFlucEfficiencyTask, processEfficiency, "Process efficiency task", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JFlucEfficiencyTask>(cfgc)};
}
