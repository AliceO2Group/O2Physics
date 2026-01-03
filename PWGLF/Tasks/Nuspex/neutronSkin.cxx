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
/// \file neutronSkin.cxx
/// \brief Task to study neutron skin effects with comprehensive efficiency analysis
/// \author DongJo Kim, Jasper Parkkila, Bong-Hwi Lim (djkim@cern.ch, jparkkil@cern.ch, bong-hwi.lim@cern.ch)
/// \since March 2024

// ============================================================================
// INCLUDES
// ============================================================================
#include "PWGLF/Utils/collisionCuts.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <string>
#include <vector>

// ============================================================================
// USING STATEMENTS
// ============================================================================
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::rctsel;
using namespace o2::constants::physics;

// ============================================================================
// TASK STRUCT
// ============================================================================
struct NeutronSkinTask {
  // ==========================================================================
  // SERVICES & CONSTANTS
  // ==========================================================================
  Service<o2::framework::O2DatabasePDG> pdg;

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

  static constexpr double ChargeThreshold = 3.0; // PDG charge units: 3 = |e|

  // ==========================================================================
  // THNSPARSE AXIS ENUMERATION
  // ==========================================================================
  enum class PtSparseAxis : int {
    kPt = 0,
    kCentrality = 1,
    kCharge = 2,
    kEtaSign = 3,
    kMcType = 4,
    kPhi = 5,
    kNAxes = 6
  };

  // ==========================================================================
  // HELPER FUNCTIONS: INDEX CONVERSION
  // ==========================================================================
  inline int getChargeIndex(int sign) const
  {
    return (sign > 0) ? 1 : 0; // 0=negative, 1=positive
  }

  inline int getEtaSignIndex(double eta) const
  {
    return (eta > 0) ? 1 : 0; // 0=negative eta, 1=positive eta
  }

  inline int getMcTypeIndex(bool isGenerated) const
  {
    return isGenerated ? 1 : 0; // 0=Reconstructed, 1=Generated
  }

  // ==========================================================================
  // CONFIGURABLES
  // ==========================================================================
  ConfigurableAxis axisPt{"axisPt", std::vector<double>(PttJacek.begin(), PttJacek.end()), "pT axis"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "multiplicity / centrality axis"};
  ConfigurableAxis axisPhi{"axisPhi", {40, 0, o2::constants::math::TwoPI}, "phi axis (rad)"};

  // Event cuts
  Configurable<int> cfgAcceptSplitCollisions{"cfgAcceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};

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
    Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", false, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", true, "Evt sel: apply pileup rejection"};
    Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtCollInTimeRangeStandard{"cfgEvtCollInTimeRangeStandard", false, "Evt sel: apply NoCollInTimeRangeStandard"};
    Configurable<bool> cfgEvtRun2AliEventCuts{"cfgEvtRun2AliEventCuts", true, "Evt sel: apply Run2 Ali event cuts"};
    Configurable<bool> cfgEvtRun2INELgtZERO{"cfgEvtRun2INELgtZERO", false, "Evt sel: apply Run2 INEL>0 event cuts"};
    Configurable<bool> cfgEvtUseRCTFlagChecker{"cfgEvtUseRCTFlagChecker", false, "Evt sel: use RCT flag checker"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
    Configurable<bool> cfgEvtGoodITSLayersAll{"cfgEvtGoodITSLayersAll", true, "Evt sel: cut time intervals with dead ITS staves"};
  } EventCuts;

  // Track selections
  struct : ConfigurableGroup {
    Configurable<float> cfgMinPt{"cfgMinPt", 0.6, "Track minimum pt cut"};
    Configurable<float> cfgMaxPt{"cfgMaxPt", 300.0f, "Maximum transverse momentum"};
    Configurable<float> cfgEtaMin{"cfgEtaMin", -1.0f, "Minimum pseudorapidity"};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 1.0f, "Maximum pseudorapidity"};
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"};
    Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", true, "Global track selection"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};
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
    Configurable<float> cfgMaxbDCArToPVcut{"cfgMaxbDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
    Configurable<float> cfgMaxbDCAzToPVcut{"cfgMaxbDCAzToPVcut", 1.0, "Track DCAz cut to PV Maximum"};
    Configurable<float> cfgPIDnSigmaCut{"cfgPIDnSigmaCut", 3.0, "PID nSigma cut"};
  } TrackCuts;

  Configurable<bool> applyMCStudy{"applyMCStudy", false, "Apply MC study"};
  Configurable<bool> debugMode{"debugMode", false, "Debug mode"};

  // ==========================================================================
  // HELPER OBJECTS
  // ==========================================================================
  o2::analysis::CollisonCuts colCuts;
  RCTFlagsChecker rctChecker;

  // ==========================================================================
  // DATA TYPE DEFINITIONS
  // ==========================================================================
  // Run3 types
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksCov>;
  using TrackCandidatesPID = soa::Join<TrackCandidates, aod::pidTPCPi>;

  // Run2 types
  using CollisionRun2Candidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  // MC types
  using MCCollisionCandidates = soa::Join<CollisionCandidates, aod::McCollisionLabel>;
  using MCTrackCandidates = soa::Join<TrackCandidates, aod::McTrackLabel>;
  using MCTrackCandidatesPID = soa::Join<MCTrackCandidates, aod::pidTPCPi>;
  using MCRun2CollisionCandidates = soa::Join<CollisionRun2Candidates, aod::McCollisionLabel>;

  // Filter declarations
  Filter cfCollisionFilter = (nabs(aod::collision::posZ) < EventCuts.cfgEvtZvtx);
  Filter cfTrackFilter = (aod::track::pt >= TrackCuts.cfgMinPt) &&
                         (aod::track::pt <= TrackCuts.cfgMaxPt) &&
                         (aod::track::eta >= TrackCuts.cfgEtaMin) &&
                         (aod::track::eta <= TrackCuts.cfgEtaMax);

  // Preslices
  Preslice<TrackCandidates> perCollision = aod::track::collisionId;
  Preslice<TrackCandidatesPID> perCollisionPID = aod::track::collisionId;

  // Histogram Registry
  HistogramRegistry registry{
    "registry",
    {{"hZVertexMC", "MC Z vertex distribution;Z vertex (cm);Centrality (%)", {HistType::kTH2F, {{200, -20, 20}, {axisMultiplicity}}}},
     {"hZVertexReco", "Reconstructed Z vertex distribution;Z vertex (cm);Centrality (%)", {HistType::kTH2F, {{200, -20, 20}, {axisMultiplicity}}}},
     {"hZVertexCorrelation", "Z vertex correlation;MC Z vertex (cm);Reco Z vertex (cm)", {HistType::kTH2F, {{200, -20, 20}, {200, -20, 20}}}}}};

  // THnSparse for pT analysis
  std::shared_ptr<THnSparse> hPtSparse; // 5D: pT × centrality × charge × eta_sign × mc_type

  // ==========================================================================
  // INITIALIZATION
  // ==========================================================================
  void init(InitContext const&)
  {
    if (debugMode) {
      LOGF(info, "Initializing NeutronSkinTask");
    }

    // Initialize collision cuts
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
    colCuts.setApplyGoodITSLayersAll(EventCuts.cfgEvtGoodITSLayersAll);
    colCuts.setApplyRun2AliEventCuts(EventCuts.cfgEvtRun2AliEventCuts);
    colCuts.setApplyRun2INELgtZERO(EventCuts.cfgEvtRun2INELgtZERO);
    colCuts.printCuts();

    // Initialize RCT checker
    rctChecker.init(EventCuts.cfgEvtRCTFlagCheckerLabel, EventCuts.cfgEvtRCTFlagCheckerZDCCheck, EventCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    // ==========================================================================
    // Initialize 6D THnSparse for pT analysis (with phi dimension)
    // ==========================================================================

    AxisSpec axisPtSparse{axisPt, "p_{T} (GeV/c)"};
    AxisSpec axisCentSparse{axisMultiplicity, "Centrality (%)"};
    AxisSpec axisChargeSparse{2, -0.5, 1.5, "Charge (0=neg, 1=pos)"};
    AxisSpec axisEtaSignSparse{2, -0.5, 1.5, "#eta sign (0=neg, 1=pos)"};
    AxisSpec axisMcTypeSparse{2, -0.5, 1.5, "Type (0=Rec, 1=Gen)"};
    AxisSpec axisPhiSparse{axisPhi, "#phi (rad)"};

    if (doprocessMC || doprocessMCRun2 || doprocessMCPID) {
      // 6D THnSparse for MC (includes mc_type and phi dimensions)
      hPtSparse = registry.add<THnSparse>(
        "hPtSparse",
        "6D p_{T} analysis;p_{T} (GeV/c);Centrality (%);Charge;#eta sign;MC Type;#phi (rad)",
        HistType::kTHnSparseD,
        {axisPtSparse, axisCentSparse, axisChargeSparse, axisEtaSignSparse, axisMcTypeSparse, axisPhiSparse});

      if (debugMode) {
        LOGF(info, "Created 6D THnSparse:");
        LOGF(info, "  pT: %d bins", hPtSparse->GetAxis(0)->GetNbins());
        LOGF(info, "  centrality: %d bins", hPtSparse->GetAxis(1)->GetNbins());
        LOGF(info, "  charge: %d bins", hPtSparse->GetAxis(2)->GetNbins());
        LOGF(info, "  eta_sign: %d bins", hPtSparse->GetAxis(3)->GetNbins());
        LOGF(info, "  mc_type: %d bins", hPtSparse->GetAxis(4)->GetNbins());
        LOGF(info, "  phi: %d bins", hPtSparse->GetAxis(5)->GetNbins());
      }
    } else {
      // 5D THnSparse for data (no mc_type dimension, but has phi)
      hPtSparse = registry.add<THnSparse>(
        "hPtSparse",
        "5D p_{T} analysis;p_{T} (GeV/c);Centrality (%);Charge;#eta sign;#phi (rad)",
        HistType::kTHnSparseD,
        {axisPtSparse, axisCentSparse, axisChargeSparse, axisEtaSignSparse, axisPhiSparse});
    }

    // Add basic eta histograms separately
    if (doprocessMC || doprocessMCRun2 || doprocessMCPID) {
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
    if (applyMCStudy && (doprocessMC || doprocessMCRun2 || doprocessMCPID)) {
      registry.add("hChargeSignMismatch", "Charge-Sign mismatch cases", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMismatchPos", "MC charge + but track sign -", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMismatchNeg", "MC charge - but track sign +", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMatch", "Charge-Sign match cases", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMatchPos", "MC charge + and track sign +", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
      registry.add("hChargeSignMatchNeg", "MC charge - and track sign -", {HistType::kTH2D, {AxisSpec(axisPt), AxisSpec(axisMultiplicity)}});
    }
  }

  // ==========================================================================
  // HELPER FUNCTIONS: PARTICLE PHYSICS
  // ==========================================================================

  /// Get the charge of a particle from PDG database
  /// \param particle MC particle
  /// \return Charge in PDG units (10.f if particle not found)
  template <typename ParticleType>
  double getCharge(ParticleType const& particle)
  {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle) {
      return 10.f;
    }
    return pdgParticle->Charge();
  }

  /// Check if PDG code corresponds to a charged particle
  bool isChargedParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= ChargeThreshold;
  }

  // ==========================================================================
  // HELPER FUNCTIONS: TRACK SELECTION
  // ==========================================================================

  /// Apply comprehensive track selection cuts
  /// \param track Reconstructed track
  /// \return true if track passes all cuts
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

  // ==========================================================================
  // HELPER FUNCTIONS: EVENT SELECTION
  // ==========================================================================

  /// Common event selection function for Run3
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

  /// Common event selection function for Run2 data
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

  // ==========================================================================
  // HELPER FUNCTIONS: HISTOGRAM FILLING
  // ==========================================================================

  /// Common histogram filling function for tracks
  template <typename TrackType>
  void fillTrackHistograms(const TrackType& track, float centrality)
  {
    // Fill THnSparse (6D for MC, 5D for data)
    fillPtSparse(track.pt(), centrality, track.sign(), track.eta(), track.phi(), false);

    // Fill eta histogram separately (continuous eta values)
    registry.fill(HIST("hEtaRec"), track.eta(), centrality);

    // pT uncertainty histograms (unchanged)
    auto ptUncertaintySigma = track.c1Pt21Pt2() * track.pt() * track.pt();
    auto ptUncertainty = std::sqrt(ptUncertaintySigma);
    registry.fill(HIST("hPtUncertainty"), track.pt(), centrality, ptUncertainty);

    if (track.sign() > 0) {
      registry.fill(HIST("hPtUncertaintyPos"), track.pt(), centrality, ptUncertainty);
    } else if (track.sign() < 0) {
      registry.fill(HIST("hPtUncertaintyNeg"), track.pt(), centrality, ptUncertainty);
    }
  }

  /// Common histogram filling function for MC particles
  template <typename ParticleType>
  void fillMCParticleHistograms(const ParticleType& particle, float centrality)
  {
    // Get particle charge
    auto charge = getCharge(particle);
    int chargeSign = (charge > 0) ? 1 : -1;

    // Fill THnSparse (6D for MC)
    fillPtSparse(particle.pt(), centrality, chargeSign, particle.eta(), particle.phi(), true);

    // Fill eta histogram separately (continuous eta values)
    registry.fill(HIST("hEtaGen"), particle.eta(), centrality);
  }

  // ==========================================================================
  // HELPER FUNCTIONS: THNSPARSE FILLING
  // ==========================================================================

  /// Fill the THnSparse for pT analysis (6D for MC, 5D for data)
  /// \param pt Transverse momentum (GeV/c)
  /// \param centrality Event centrality (%)
  /// \param charge Track charge sign (>0 for positive, <0 for negative)
  /// \param eta Pseudorapidity
  /// \param phi Azimuthal angle (rad)
  /// \param isGenerated true for MC generated, false for reconstructed
  void fillPtSparse(double pt, float centrality, int charge, double eta, double phi, bool isGenerated = false)
  {
    if (!hPtSparse) {
      LOGF(error, "hPtSparse not initialized!");
      return;
    }

    int chargeIdx = getChargeIndex(charge);
    int etaSignIdx = getEtaSignIndex(eta);

    if (doprocessMC || doprocessMCRun2 || doprocessMCPID) {
      // 6D filling for MC processing (includes mc_type and phi)
      int mcTypeIdx = getMcTypeIndex(isGenerated);
      double fillArray[6] = {
        pt,
        centrality,
        static_cast<double>(chargeIdx),
        static_cast<double>(etaSignIdx),
        static_cast<double>(mcTypeIdx),
        phi};
      hPtSparse->Fill(fillArray);
    } else {
      // 5D filling for data processing (no mc_type, but has phi)
      double fillArray[5] = {
        pt,
        centrality,
        static_cast<double>(chargeIdx),
        static_cast<double>(etaSignIdx),
        phi};
      hPtSparse->Fill(fillArray);
    }
  }

  // ==========================================================================
  // PROCESS METHODS
  // ==========================================================================

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

  void processMC(aod::McCollisions::iterator const& mcCollision,
                 soa::SmallGroups<MCCollisionCandidates> const& collisions,
                 soa::Filtered<MCTrackCandidates> const& mcTracks,
                 aod::McParticles const& mcParticles)
  {
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

  void processMCPID(aod::McCollisions::iterator const& mcCollision,
                    soa::SmallGroups<MCCollisionCandidates> const& collisions,
                    soa::Filtered<MCTrackCandidatesPID> const& mcTracks,
                    aod::McParticles const& mcParticles)
  {
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

  PROCESS_SWITCH(NeutronSkinTask, processMC, "Process MC only", false);
  PROCESS_SWITCH(NeutronSkinTask, processMCPID, "Process MC with PID only", false);
  PROCESS_SWITCH(NeutronSkinTask, processMCRun2, "Process Run2 MC only", false);
  PROCESS_SWITCH(NeutronSkinTask, processData, "Process data only", false);
  PROCESS_SWITCH(NeutronSkinTask, processDataPID, "Process data with PID only", false);
  PROCESS_SWITCH(NeutronSkinTask, processDataRun2, "Process Run2 data only", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NeutronSkinTask>(cfgc)};
}
