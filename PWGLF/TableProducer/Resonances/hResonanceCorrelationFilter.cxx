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
/// \brief This task pre-filters tracks to do h-resonance (phi, K*0)
///        correlations with an analysis task.
///
/// \file hResonanceCorrelationFilter.cxx
/// \author Hirak Kumar Koley (hirak.koley@cern.ch)

#include "PWGLF/DataModel/LFHResonanceCorrelationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TList.h>
#include <TPDGCode.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::constants::math;
using namespace o2::aod::rctsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

using LorentzVectorPtEtaPhiMass = ROOT::Math::PtEtaPhiMVector;

enum PIDCutType {
  SquareType = 1,
  CircularType,
};

#define BIT_SET(var, nbit) ((var) |= (1 << (nbit)))
#define BIT_CHECK(var, nbit) ((var) & (1 << (nbit)))

struct HResonanceCorrelationFilter {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  RCTFlagsChecker rctChecker;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};

  // Operational
  Configurable<std::string> zorroMask{"zorroMask", "", "zorro trigger class to select on (empty: none)"};

  Configurable<float> cfgTPCNsigmaKaon{"cfgTPCNsigmaKaon", 3.0f, "TPC Kaon PID"};
  Configurable<float> cfgTOFNsigmaKaon{"cfgTOFNsigmaKaon", 3.0f, "TOF Kaon PID"};

  Configurable<float> cfgMinMass{"cfgMinMass", 1.005f, "Minimum KK mass"};
  Configurable<float> cfgMaxMass{"cfgMaxMass", 1.035f, "Maximum KK mass"};

  Configurable<float> cfgMinMassKstar{"cfgMinMassKstar", 0.796f, "Minimum K#pi mass"};
  Configurable<float> cfgMaxMassKstar{"cfgMaxMassKstar", 0.996f, "Maximum K#pi mass"};

  Configurable<float> cfgRapidity{"cfgRapidity", 0.5f, "Rapidity cut"};

  // used for event selections in Pb-Pb
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low cut on TPC occupancy"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections";
    // event filtering
    Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
    Configurable<bool> selectINELgtZERO{"selectINELgtZERO", true, "select INEL>0 events"};
    Configurable<bool> requireAllGoodITSLayers{"requireAllGoodITSLayers", false, " require that in the event all ITS are good"};
    Configurable<bool> requireGoodTriggerTVX{"requireGoodTriggerTVX", false, " require acceptable FT0C-FT0A time difference"};
    Configurable<bool> requireGoodZvtxFT0vsPV{"requireGoodZvtxFT0vsPV", false, " require small difference between z-vertex from PV and from FT0"};
    Configurable<float> minCentPercent{"minCentPercent", 0, "minimum centrality percentage"};
    Configurable<float> maxCentPercent{"maxCentPercent", 100, "maximum centrality percentage"};
  } eventSelections;

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
    std::string prefix = "generalSelections";

    // Associated particle selections in phase space
    Configurable<float> assocEtaMin{"assocEtaMin", -0.8, "triggeretamin"};
    Configurable<float> assocEtaMax{"assocEtaMax", 0.8, "triggeretamax"};
    Configurable<float> assocPtCutMin{"assocPtCutMin", 0.2, "assocptmin"};
    Configurable<float> assocPtCutMax{"assocPtCutMax", 10, "assocptmax"};

    // Trigger particle selections in phase space
    Configurable<float> triggerEtaMin{"triggerEtaMin", -0.8, "triggeretamin"};
    Configurable<float> triggerEtaMax{"triggerEtaMax", 0.8, "triggeretamax"};
    Configurable<float> triggerPtCutMin{"triggerPtCutMin", 3, "triggerptmin"};
    Configurable<float> triggerPtCutMax{"triggerPtCutMax", 20, "triggerptmax"};
  } generalSelections;

  struct : ConfigurableGroup {
    std::string prefix = "trackSelections";
    // Track quality
    Configurable<int> minTPCNCrossedRows{"minTPCNCrossedRows", 70, "Minimum TPC crossed rows"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<bool> assocRequireITS{"assocRequireITS", true, "require ITS signal in assoc tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};
    Configurable<bool> requireClusterInITS{"requireClusterInITS", false, "require cluster in ITS for phi daughter tracks"};
    Configurable<int> minITSClustersForDaughterTracks{"minITSClustersForDaughterTracks", 1, "Minimum number of ITS clusters for phi daughter tracks"};

    // Associated pion identification
    Configurable<float> pionMinBayesProb{"pionMinBayesProb", 0.95, "minimal Bayesian probability for pion ID"};
    Configurable<float> assocPionNSigmaTPCFOF{"assocPionNSigmaTPCFOF", 3, "minimal n sigma in TOF and TPC for Pion ID"};
    Configurable<float> rejectSigma{"rejectSigma", 1, "n sigma for rejecting pion candidates"};

    // primary particle DCAxy selections
    // formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};
  } trackSelections;

  struct : ConfigurableGroup {
    /// PID Selections
    Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 4.0f, "pidnSigma Cut for pre-selection of tracks"};
    Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};                       // By pass TOF PID selection
    Configurable<int> cPIDcutType{"cPIDcutType", 2, "cPIDcutType = 1 for square cut, 2 for circular cut"}; // By pass TOF PID selection
    Configurable<bool> ispTdepPID{"ispTdepPID", false, "enable pT dependent PID"};

    // Kaon
    Configurable<std::vector<float>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {0.5f}, "pT intervals for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
    Configurable<std::vector<float>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.0f}, "pT intervals for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
    Configurable<std::vector<float>> kaonTPCTOFCombinedpTintv{"kaonTPCTOFCombinedpTintv", {999.0f}, "pT intervals for Kaon TPC-TOF PID cuts"};
    Configurable<std::vector<float>> kaonTPCTOFCombinedPIDcuts{"kaonTPCTOFCombinedPIDcuts", {2}, "nSigma list for Kaon TPC-TOF PID cuts"};

    // Pion (K*0 daughter)
    Configurable<std::vector<float>> pionTPCPIDpTintv{"pionTPCPIDpTintv", {0.5f}, "pT intervals for Pion TPC PID cuts"};
    Configurable<std::vector<float>> pionTPCPIDcuts{"pionTPCPIDcuts", {2}, "nSigma list for Pion TPC PID cuts"};
    Configurable<std::vector<float>> pionTOFPIDpTintv{"pionTOFPIDpTintv", {999.0f}, "pT intervals for Pion TOF PID cuts"};
    Configurable<std::vector<float>> pionTOFPIDcuts{"pionTOFPIDcuts", {2}, "nSigma list for Pion TOF PID cuts"};
    Configurable<std::vector<float>> pionTPCTOFCombinedpTintv{"pionTPCTOFCombinedpTintv", {999.0f}, "pT intervals for Pion TPC-TOF PID cuts"};
    Configurable<std::vector<float>> pionTPCTOFCombinedPIDcuts{"pionTPCTOFCombinedPIDcuts", {2}, "nSigma list for Pion TPC-TOF PID cuts"};
  } configPID;

  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> parameterCCDBPath{"parameterCCDBPath", "Users/k/kcui/LHC25b4a/parameter", "Path of the mean and sigma"};

  // must include windows for background and peak
  Configurable<float> maxMassNSigma{"maxMassNSigma", 12.0f, "max mass region to be considered for further analysis"};

  // For extracting strangeness mass QA plots
  struct : ConfigurableGroup {
    ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    ConfigurableAxis axisPhiMass{"axisPhiMass", {200, 0.99f, 1.08f}, "M(K^{+}K^{-})"};
    ConfigurableAxis axisKstarMass{"axisKstarMass", {200, 0.75f, 1.05f}, "M(K#pi)"};
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Centrality percentile bins"};
  } axesConfigurations;

  // QA
  Configurable<bool> doTrueSelectionInMass{"doTrueSelectionInMass", false, "Fill mass histograms only with true primary Particles for MC"};
  // Do declarative selections for DCAs, if possible
  Filter preFilterTracks = nabs(aod::track::dcaXY) < trackSelections.dcaXYconstant + trackSelections.dcaXYpTdep * nabs(aod::track::signed1Pt);

  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  using FullTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA>;
  using DauTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA, aod::McTrackLabels>;
  // using IDTracks= soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::TOFSignal>; // prepared for Bayesian PID
  using IDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::TOFSignal, aod::TracksDCA>;
  using IDTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::TOFSignal, aod::TracksDCA, aod::McTrackLabels>;

  Produces<aod::TriggerTracks> triggerTrack;
  Produces<aod::TriggerTrackExtras> triggerTrackExtra;
  Produces<aod::AssocPhis> assocPhis;
  Produces<aod::AssocKstars> assocKstars;
  Produces<aod::AssocHadrons> assocHadrons;
  Produces<aod::AssocPID> assocPID;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::FullTracks, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>>;

  // for MC reco
  using MCEventCandidates = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using MCTrackCandidates = soa::Filtered<soa::Join<TrackCandidates, aod::McTrackLabels>>;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  int mRunNumber = -1;

  struct TriggCandidate {
    float pt = 0.f;
    int collisionId = -1;
    int trackId = -1;
    bool isPhysicalPrimary = false;
    float origPt = 0.f;
  };
  TriggCandidate thisTrigg;

  std::vector<TriggCandidate> triggerCandidates;

  void init(InitContext const&)
  {
    rctChecker.init(configEvents.cfgEvtRCTFlagCheckerLabel, configEvents.cfgEvtRCTFlagCheckerZDCCheck, configEvents.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

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

    zorroSummary.setObject(zorro.getZorroSummary());
    mRunNumber = -1;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), zorroMask.value);
    zorro.populateHistRegistry(histos, bc.runNumber());

    mRunNumber = bc.runNumber();
  }

  template <typename Coll>
  bool isSelected(const Coll& collision, bool fillHist = true)
  {
    auto applyCut = [&](bool enabled, bool condition, int bin) {
      if (!enabled) {
        return true;
      }
      if (!condition) {
        return false;
      }
      if (fillHist) {
        histos.fill(HIST("CollCutCounts"), bin);
      }
      return true;
    };

    if (fillHist) {
      histos.fill(HIST("CollCutCounts"), 0);
    }

    if (!applyCut(true, std::abs(collision.posZ()) <= configEvents.cfgEvtZvtx, 1)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtTriggerTVXSel,
                  collision.selection_bit(aod::evsel::kIsTriggerTVX), 2)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtNoTFBorderCut,
                  collision.selection_bit(aod::evsel::kNoTimeFrameBorder), 3)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtNoITSROFrameBorderCut,
                  collision.selection_bit(aod::evsel::kNoITSROFrameBorder), 4)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtIsRCTFlagpassed, rctChecker(collision), 5)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtSel8, collision.sel8(), 6)) {
      return false;
    }

    if (!applyCut(configEvents.cfgEvtIsINELgt0, collision.isInelGt0(), 7)) {
      return false;
    }

    if (fillHist) {
      histos.fill(HIST("CollCutCounts"), 8);
    }

    return true;
  }

  // this function allows for all event selections to be done in a modular way
  template <typename TCollision>
  bool isCollisionSelected(TCollision const& collision)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return false;
    }
    if (std::abs(collision.posZ()) > eventSelections.zVertexCut) {
      return false;
    }
    if (collision.centFT0M() > eventSelections.maxCentPercent || collision.centFT0M() < eventSelections.minCentPercent) {
      return false;
    }
    if (!collision.isInelGt0() && eventSelections.selectINELgtZERO) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodITSLayersAll) && eventSelections.requireAllGoodITSLayers) {
      return false;
    }
    if (zorroMask.value != "") {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return false;
      }
    }
    return true;
  }

  // more event selections in Pb-Pb
  template <typename TCollision>
  bool isCollisionSelectedPbPb(TCollision collision)
  {
    if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) && eventSelections.requireGoodTriggerTVX) /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) && eventSelections.requireAllGoodITSLayers) // cut time intervals with dead ITS staves
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) && eventSelections.requireGoodZvtxFT0vsPV) // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      return false;
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh) /* Below min occupancy and Above max occupancy*/
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) // reject collisions close to Time Frame borders
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) // reject events affected by the ITS ROF border
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      return false;
    return true;
  }

  template <typename TrackType>
  bool trackCut(const TrackType& track)
  {
    // basic track cuts
    if (configTracks.cDCAr7SigCut && std::abs(track.dcaXY()) > (0.004f + 0.013f / (track.pt()))) // 7 - Sigma cut
    {
      return false;
    }
    if (configTracks.cTPCNClsFound && (track.tpcNClsFound() < configTracks.cMinTPCNClsFound)) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < configTracks.cfgMinCrossedRows) {
      return false;
    }
    if (configTracks.cfgHasTOF && !track.hasTOF()) {
      return false;
    }
    if (configTracks.cfgPrimaryTrack && !track.isPrimaryTrack()) {
      return false;
    }
    if (configTracks.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (configTracks.cfgPVContributor && !track.isPVContributor()) {
      return false;
    }
    if (configTracks.cfgGlobalTrack && !track.isGlobalTrack()) {
      return false;
    }

    return true;
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

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    auto vKaonTPCPIDcuts = configPID.kaonTPCPIDcuts.value;
    auto vKaonTPCTOFCombinedPIDcuts = configPID.kaonTPCTOFCombinedPIDcuts.value;

    if (!configPID.cByPassTOF && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (vKaonTPCTOFCombinedPIDcuts[0] * vKaonTPCTOFCombinedPIDcuts[0])) {
      return true;
    }
    if (!configPID.cByPassTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < vKaonTPCPIDcuts[0]) {
      return true;
    }
    if (configPID.cByPassTOF && std::abs(candidate.tpcNSigmaKa()) < vKaonTPCPIDcuts[0]) {
      return true;
    }
    return false;
  }

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

    // TPC PID interval-based check
    for (size_t i = 0; i < vPionTPCPIDpTintv.size() - 1; ++i) {
      if (pt > vPionTPCPIDpTintv[i] && pt < vPionTPCPIDpTintv[i + 1]) {
        if (std::abs(tpcNsigmaPi) < vPionTPCPIDcuts[i]) {
          tpcPIDPassed = true;
          break;
        }
      }
    }

    // TOF bypass option
    if (configPID.cByPassTOF) {
      return std::abs(tpcNsigmaPi) < vPionTPCPIDcuts.back();
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
        for (size_t i = 0; i < vPionTOFPIDpTintv.size(); ++i) {
          if (pt < vPionTOFPIDpTintv[i]) {
            if (std::abs(tofNsigmaPi) < vPionTOFPIDcuts[i] &&
                std::abs(tpcNsigmaPi) < vPionTPCPIDcuts.back()) {
              return true;
            }
          }
        }
      } else if (configPID.cPIDcutType == CircularType) {
        // Circular cut
        for (size_t i = 0; i < vPionTPCTOFCombinedpTintv.size(); ++i) {
          if (pt < vPionTPCTOFCombinedpTintv[i]) {
            float combinedSigma2 = tpcNsigmaPi * tpcNsigmaPi +
                                   tofNsigmaPi * tofNsigmaPi;
            if (combinedSigma2 < vPionTPCTOFCombinedPIDcuts[i] * vPionTPCTOFCombinedPIDcuts[i]) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    auto vPionTPCPIDcuts = configPID.pionTPCPIDcuts.value;
    auto vPionTPCTOFCombinedPIDcuts = configPID.pionTPCTOFCombinedPIDcuts.value;

    if (!configPID.cByPassTOF && candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (vPionTPCTOFCombinedPIDcuts[0] * vPionTPCTOFCombinedPIDcuts[0])) {
      return true;
    }
    if (!configPID.cByPassTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < vPionTPCPIDcuts[0]) {
      return true;
    }
    if (configPID.cByPassTOF && std::abs(candidate.tpcNSigmaPi()) < vPionTPCPIDcuts[0]) {
      return true;
    }
    return false;
  }

  // reco-level trigger quality checks (N.B.: DCA is filtered, not selected)
  template <class TTrack>
  bool isValidTrigger(TTrack track)
  {
    if (track.eta() > generalSelections.triggerEtaMax || track.eta() < generalSelections.triggerEtaMin) {
      return false;
    }
    // if (track.sign()= 1 ) {continue;}
    if (track.pt() > generalSelections.triggerPtCutMax || track.pt() < generalSelections.triggerPtCutMin) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < trackSelections.minTPCNCrossedRows) {
      return false; // crossed rows
    }
    if (!track.hasITS() && trackSelections.triggerRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > trackSelections.triggerMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(BIT_CHECK(track.itsClusterMap(), 0)) && trackSelections.triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    return true;
  }

  template <class TTrack>
  bool isValidAssocTrack(TTrack assoc)
  {
    if (assoc.eta() > generalSelections.assocEtaMax || assoc.eta() < generalSelections.assocEtaMin) {
      return false;
    }
    if (assoc.pt() > generalSelections.assocPtCutMax || assoc.pt() < generalSelections.assocPtCutMin) {
      return false;
    }
    if (assoc.tpcNClsCrossedRows() < trackSelections.minTPCNCrossedRows) {
      return false; // crossed rows
    }
    if (!assoc.hasITS() && trackSelections.assocRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }

    // do this only if information is available
    float nSigmaTPCTOF[8] = {-10, -10, -10, -10, -10, -10, -10, -10};
    if constexpr (requires { assoc.tofSignal(); } && !requires { assoc.mcParticle(); }) {
      if (assoc.tofSignal() > 0) {
        if (std::sqrt(assoc.tofNSigmaPi() * assoc.tofNSigmaPi() + assoc.tpcNSigmaPi() * assoc.tpcNSigmaPi()) > trackSelections.assocPionNSigmaTPCFOF)
          return false;
        if (assoc.tofNSigmaPr() < trackSelections.rejectSigma)
          return false;
        if (assoc.tpcNSigmaPr() < trackSelections.rejectSigma)
          return false;
        if (assoc.tofNSigmaKa() < trackSelections.rejectSigma)
          return false;
        if (assoc.tpcNSigmaKa() < trackSelections.rejectSigma)
          return false;
        nSigmaTPCTOF[4] = assoc.tofNSigmaPi();
        nSigmaTPCTOF[5] = assoc.tofNSigmaKa();
        nSigmaTPCTOF[6] = assoc.tofNSigmaPr();
        nSigmaTPCTOF[7] = assoc.tofNSigmaEl();
      } else {
        if (assoc.tpcNSigmaPi() > trackSelections.assocPionNSigmaTPCFOF)
          return false;
        if (assoc.tpcNSigmaPr() < trackSelections.rejectSigma)
          return false;
        if (assoc.tpcNSigmaKa() < trackSelections.rejectSigma)
          return false;
      }
      nSigmaTPCTOF[0] = assoc.tpcNSigmaPi();
      nSigmaTPCTOF[1] = assoc.tpcNSigmaKa();
      nSigmaTPCTOF[2] = assoc.tpcNSigmaPr();
      nSigmaTPCTOF[3] = assoc.tpcNSigmaEl();
    }

    bool physicalPrimary = false;
    float origPt = -1;
    float code = -9999;
    if constexpr (requires { assoc.mcParticle(); }) {
      if (assoc.has_mcParticle()) {
        auto mcParticle = assoc.mcParticle();
        physicalPrimary = mcParticle.isPhysicalPrimary();
        origPt = mcParticle.pt();
        code = mcParticle.pdgCode();
      }
    }

    assocHadrons(
      assoc.collisionId(),
      physicalPrimary,
      assoc.globalIndex(),
      origPt,
      code);
    assocPID(
      nSigmaTPCTOF[0],
      nSigmaTPCTOF[1],
      nSigmaTPCTOF[2],
      nSigmaTPCTOF[3],
      nSigmaTPCTOF[4],
      nSigmaTPCTOF[5],
      nSigmaTPCTOF[6],
      nSigmaTPCTOF[7]);
    return true;
  }

  // for real data processing
  void processTriggers(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, soa::Filtered<FullTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    triggerCandidates.clear();
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    double leadingPt = -1.;
    int leadingId = -1;
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      thisTrigg.pt = track.pt();
      thisTrigg.trackId = track.globalIndex();
      thisTrigg.collisionId = track.collisionId();
      thisTrigg.isPhysicalPrimary = false; // if you decide to check real data for primaries, you'll have a hard time
      thisTrigg.origPt = 0;
      triggerCandidates.push_back(thisTrigg);
      if (track.pt() > leadingPt) {
        leadingPt = track.pt();
        leadingId = track.globalIndex();
      }
    }
    for (auto const& TriggCandidate : triggerCandidates) {
      bool isLeading = (leadingId == TriggCandidate.trackId);
      triggerTrack(
        TriggCandidate.collisionId,
        TriggCandidate.isPhysicalPrimary,
        TriggCandidate.trackId,
        TriggCandidate.origPt,
        isLeading);
      triggerTrackExtra(1);
    }
  }

  // for MC processing
  void processTriggersMC(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, soa::Filtered<FullTracksMC> const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    triggerCandidates.clear();
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    double leadingPt = -1.;
    int leadingId = -1;
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      thisTrigg.pt = track.pt();
      thisTrigg.trackId = track.globalIndex();
      thisTrigg.collisionId = track.collisionId();
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        thisTrigg.isPhysicalPrimary = mcParticle.isPhysicalPrimary();
        thisTrigg.origPt = mcParticle.pt();
      }
      triggerCandidates.push_back(thisTrigg);
      if (track.pt() > leadingPt) {
        leadingPt = track.pt();
        leadingId = track.globalIndex();
      }
    }

    for (auto const& TriggCandidate : triggerCandidates) {
      bool isLeading = (leadingId == TriggCandidate.trackId);
      triggerTrack(
        TriggCandidate.collisionId,
        TriggCandidate.isPhysicalPrimary,
        TriggCandidate.trackId,
        TriggCandidate.origPt,
        isLeading);
      triggerTrackExtra(1);
    }
  }

  void processAssocPions(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > eventSelections.zVertexCut) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }

  void processAssocPionsMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracksMC> const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > eventSelections.zVertexCut) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }

  void processAssocHadrons(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<FullTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > eventSelections.zVertexCut) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }
  void processAssocHadronsMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<FullTracksMC> const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > eventSelections.zVertexCut) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }

  void processPhis(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks)
  {
    if (!isSelected(collision)) {
      return;
    }

    LorentzVectorPtEtaPhiMass kplus;
    LorentzVectorPtEtaPhiMass kminus;
    LorentzVectorPtEtaPhiMass phi;

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      // same collision
      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      // quality cuts
      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      // unlike-sign only
      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      // kaon PID
      if (configPID.ispTdepPID) {
        if (!ptDependentPidKaon(trk1) ||
            !ptDependentPidKaon(trk2)) {
          continue;
        }
      } else {
        if (!selectionPID(trk1) ||
            !selectionPID(trk2)) {
          continue;
        }
      }

      // assign K+ and K-
      if (trk1.sign() > 0) {
        kplus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                           o2::constants::physics::MassKPlus);
      } else {
        kplus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                           o2::constants::physics::MassKPlus);
      }

      phi = kplus + kminus;

      float invMass = phi.M();

      // rapidity cut
      if (std::abs(phi.Rapidity()) > configTracks.cfgCutRapidity) {
        continue;
      }

      // optional mass window
      if (invMass < cfgMinMass || invMass > cfgMaxMass) {
        continue;
      }

      assocPhis(
        collision.globalIndex(),
        false,
        false,
        phi.Pt(),
        phi.Eta(),
        phi.Phi(),
        invMass,
        trk1.globalIndex(),
        trk2.globalIndex());
    }
  }

  void processPhisMC(
    MCEventCandidates::iterator const& collision,
    aod::McCollisions const&,
    MCTrackCandidates const& tracks,
    aod::McParticles const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!isSelected(collision)) {
      return;
    }

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      // -----------------------------
      // Same cuts as phi1020analysis
      // -----------------------------

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      if (!selectionPID(trk1) || !selectionPID(trk2)) {
        continue;
      }

      LorentzVectorPtEtaPhiMass kplus;
      LorentzVectorPtEtaPhiMass kminus;
      LorentzVectorPtEtaPhiMass phi;

      // assign K+ and K-
      if (trk1.sign() > 0) {
        kplus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                           o2::constants::physics::MassKPlus);
      } else {
        kplus = LorentzVectorPtEtaPhiMass(trk2.pt(), trk2.eta(), trk2.phi(),
                                          o2::constants::physics::MassKPlus);
        kminus = LorentzVectorPtEtaPhiMass(trk1.pt(), trk1.eta(), trk1.phi(),
                                           o2::constants::physics::MassKPlus);
      }

      phi = kplus + kminus;

      float invMass = phi.M();

      // rapidity cut
      if (std::abs(phi.Rapidity()) > configTracks.cfgCutRapidity) {
        continue;
      }

      // optional mass window
      if (invMass < cfgMinMass || invMass > cfgMaxMass) {
        continue;
      }

      bool mcTruePhi = false;
      bool mcPhysicalPrimary = false;

      // -----------------------------
      // MC matching
      // -----------------------------

      if (trk1.has_mcParticle() && trk2.has_mcParticle()) {

        auto mc1 = trk1.mcParticle();
        auto mc2 = trk2.mcParticle();

        if (mc1.has_mothers() && mc2.has_mothers()) {

          for (auto const& mother1 : mc1.mothers_as<aod::McParticles>()) {
            for (auto const& mother2 : mc2.mothers_as<aod::McParticles>()) {

              if (mother1.globalIndex() != mother2.globalIndex()) {
                continue;
              }

              if (std::abs(mother1.pdgCode()) != 333) {
                continue;
              }

              mcTruePhi = true;
              mcPhysicalPrimary = mother1.isPhysicalPrimary();
            }
          }
        }
      }

      assocPhis(
        collision.globalIndex(),
        mcTruePhi,
        mcPhysicalPrimary,
        phi.Pt(),
        phi.Eta(),
        phi.Phi(),
        invMass,
        trk1.globalIndex(),
        trk2.globalIndex());
    }
  }

  // K*0 -> K+ pi- (and c.c.) reconstruction. Unlike phi (K+K-, symmetric daughter
  // mass), the two daughter tracks are different species, so both mass hypotheses
  // (trk1=kaon,trk2=pion) and (trk1=pion,trk2=kaon) are tried independently; either,
  // both, or neither may pass PID for a given unlike-sign pair.
  void processKstars(EventCandidates::iterator const& collision,
                     TrackCandidates const& tracks)
  {
    if (!isSelected(collision)) {
      return;
    }

    LorentzVectorPtEtaPhiMass kaon;
    LorentzVectorPtEtaPhiMass pion;
    LorentzVectorPtEtaPhiMass kstar;

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      // same collision
      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      // quality cuts
      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      // unlike-sign only
      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      for (int hypothesis = 0; hypothesis < 2; hypothesis++) {
        auto const& kaonTrack = (hypothesis == 0) ? trk1 : trk2;
        auto const& pionTrack = (hypothesis == 0) ? trk2 : trk1;

        // kaon + pion PID
        bool passKaonPID = configPID.ispTdepPID ? ptDependentPidKaon(kaonTrack) : selectionPID(kaonTrack);
        bool passPionPID = configPID.ispTdepPID ? ptDependentPidPion(pionTrack) : selectionPIDPion(pionTrack);
        if (!passKaonPID || !passPionPID) {
          continue;
        }

        kaon = LorentzVectorPtEtaPhiMass(kaonTrack.pt(), kaonTrack.eta(), kaonTrack.phi(),
                                         o2::constants::physics::MassKPlus);
        pion = LorentzVectorPtEtaPhiMass(pionTrack.pt(), pionTrack.eta(), pionTrack.phi(),
                                         o2::constants::physics::MassPiPlus);

        kstar = kaon + pion;

        float invMass = kstar.M();

        // rapidity cut
        if (std::abs(kstar.Rapidity()) > configTracks.cfgCutRapidity) {
          continue;
        }

        // mass window
        if (invMass < cfgMinMassKstar || invMass > cfgMaxMassKstar) {
          continue;
        }

        auto const& posTrack = (kaonTrack.sign() > 0) ? kaonTrack : pionTrack;
        auto const& negTrack = (kaonTrack.sign() > 0) ? pionTrack : kaonTrack;

        assocKstars(
          collision.globalIndex(),
          false,
          false,
          kstar.Pt(),
          kstar.Eta(),
          kstar.Phi(),
          invMass,
          posTrack.globalIndex(),
          negTrack.globalIndex());
      }
    }
  }

  void processKstarsMC(
    MCEventCandidates::iterator const& collision,
    aod::McCollisions const&,
    MCTrackCandidates const& tracks,
    aod::McParticles const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!isSelected(collision)) {
      return;
    }

    for (auto const& [trk1, trk2] :
         combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      // -----------------------------
      // Same cuts as processKstars
      // -----------------------------

      if (trk2.index() <= trk1.index()) {
        continue;
      }

      if (trk1.collisionId() != collision.globalIndex() ||
          trk2.collisionId() != collision.globalIndex()) {
        continue;
      }

      if (!trackCut(trk1) || !trackCut(trk2)) {
        continue;
      }

      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      for (int hypothesis = 0; hypothesis < 2; hypothesis++) {
        auto const& kaonTrack = (hypothesis == 0) ? trk1 : trk2;
        auto const& pionTrack = (hypothesis == 0) ? trk2 : trk1;

        if (!selectionPID(kaonTrack) || !selectionPIDPion(pionTrack)) {
          continue;
        }

        LorentzVectorPtEtaPhiMass kaon = LorentzVectorPtEtaPhiMass(kaonTrack.pt(), kaonTrack.eta(), kaonTrack.phi(),
                                                                   o2::constants::physics::MassKPlus);
        LorentzVectorPtEtaPhiMass pion = LorentzVectorPtEtaPhiMass(pionTrack.pt(), pionTrack.eta(), pionTrack.phi(),
                                                                   o2::constants::physics::MassPiPlus);
        LorentzVectorPtEtaPhiMass kstar = kaon + pion;

        float invMass = kstar.M();

        // rapidity cut
        if (std::abs(kstar.Rapidity()) > configTracks.cfgCutRapidity) {
          continue;
        }

        // mass window
        if (invMass < cfgMinMassKstar || invMass > cfgMaxMassKstar) {
          continue;
        }

        bool mcTrueKstar = false;
        bool mcPhysicalPrimary = false;

        // -----------------------------
        // MC matching: common mother with |pdgCode| == 313 (K*0/anti-K*0)
        // -----------------------------

        if (kaonTrack.has_mcParticle() && pionTrack.has_mcParticle()) {

          auto mcKaon = kaonTrack.mcParticle();
          auto mcPion = pionTrack.mcParticle();

          if (mcKaon.has_mothers() && mcPion.has_mothers()) {

            for (auto const& motherKaon : mcKaon.mothers_as<aod::McParticles>()) {
              for (auto const& motherPion : mcPion.mothers_as<aod::McParticles>()) {

                if (motherKaon.globalIndex() != motherPion.globalIndex()) {
                  continue;
                }

                if (std::abs(motherKaon.pdgCode()) != 313) {
                  continue;
                }

                mcTrueKstar = true;
                mcPhysicalPrimary = motherKaon.isPhysicalPrimary();
              }
            }
          }
        }

        auto const& posTrack = (kaonTrack.sign() > 0) ? kaonTrack : pionTrack;
        auto const& negTrack = (kaonTrack.sign() > 0) ? pionTrack : kaonTrack;

        assocKstars(
          collision.globalIndex(),
          mcTrueKstar,
          mcPhysicalPrimary,
          kstar.Pt(),
          kstar.Eta(),
          kstar.Phi(),
          invMass,
          posTrack.globalIndex(),
          negTrack.globalIndex());
      }
    }
  }

  PROCESS_SWITCH(HResonanceCorrelationFilter, processTriggers, "Produce trigger tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processTriggersMC, "Produce trigger tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocPions, "Produce associated Pion tables", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocPionsMC, "Produce associated Pion tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocHadrons, "Produce associated Hadron tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processAssocHadronsMC, "Produce associated Hadron tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processPhis, "Produce associated phi tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processPhisMC, "Produce associated phi tables for MC", false);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processKstars, "Produce associated K*0 tables", true);
  PROCESS_SWITCH(HResonanceCorrelationFilter, processKstarsMC, "Produce associated K*0 tables for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HResonanceCorrelationFilter>(cfgc)};
}
