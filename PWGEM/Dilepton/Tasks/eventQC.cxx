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
//
// ========================
//
// This code is for event QC for PWG-EM.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"

#include "TString.h"

#include <algorithm>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct eventQC {
  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorFV0AVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MyCollisions_Qvec = soa::Join<MyCollisions, MyQvectors>;

  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                             aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                             aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
  using MyTrack = MyTracks::iterator;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<std::string> cfg_swt_name{"cfg_swt_name", "fHighTrackMult", "desired software trigger name"};
  Configurable<bool> cfgFillEvent{"cfgFillEvent", false, "fill event histograms"};
  Configurable<bool> cfgFillTrack{"cfgFillTrack", false, "fill track histograms"};
  Configurable<bool> cfgFillPID{"cfgFillPID", false, "fill PID histograms"};
  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2, 3}, "Modulation of interest. Please keep increasing order"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
  Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
  ConfigurableAxis ConfPtBins{"ConfPtBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pT bins for output histograms"};
  Configurable<int> cfgNbinsEta{"cfgNbinsEta", 20, "number of eta bins for output histograms"};
  Configurable<int> cfgNbinsPhi{"cfgNbinsPhi", 36, "number of phi bins for output histograms"};

  ConfigurableAxis ConfFT0AMultBins{"ConfFT0AMultBins", {200, 0, 200e+3}, "FT0A multiplicity bins for output histograms"};
  ConfigurableAxis ConfFT0CMultBins{"ConfFT0CMultBins", {600, 0, 60e+3}, "FT0C multiplicity bins for output histograms"};
  ConfigurableAxis ConfFV0AMultBins{"ConfFV0AMultBins", {200, 0, 200e+3}, "FV0A multiplicity bins for output histograms"};
  ConfigurableAxis ConfTrackMultBins{"ConfTrackMultBins", {6001, -0.5, 6e+3 + 0.5}, "Track multiplicity bins for output histograms"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {110, 0, 110}, "centrality bins for output histograms"};

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -1e+10, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 1e+10, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", false, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequirekNoCollInRofStandard{"cfgRequirekNoCollInRofStandard", false, "require no other collisions in this Readout Frame with per-collision multiplicity above threshold"};
    Configurable<bool> cfgRequirekNoCollInRofStrict{"cfgRequirekNoCollInRofStrict", false, "require no other collisions in this Readout Frame"};
    Configurable<bool> cfgRequirekNoHighMultCollInPrevRof{"cfgRequirekNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
    Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "number of inactive chips on ITS layer 3 are below threshold "};
    Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "number of inactive chips on ITS layers 0-3 are below threshold "};
    Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "number of inactive chips on all ITS layers are below threshold "};
  } eventcuts;

  struct : ConfigurableGroup {
    std::string prefix = "trackcut_group";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.15, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 80, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2/NclsTOF"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.2, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.2, "max dca Z for single track in cm"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -1e+10, "min n sigma e in TPC"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +1e+10, "max n sigma e in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", 0.0, "min n sigma pi in TPC for exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", 0.0, "max n sigma pi in TPC for exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -1e+10, "min n sigma e in TOF"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +1e+10, "max n sigma e in TOF"};
    Configurable<bool> cfg_requireTOF{"cfg_requireTOF", false, "require TOF hit"};
  } trackcuts;

  // for RCT
  Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
  Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
  Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  Zorro zorro;
  std::vector<int> mTOIidx;
  uint64_t mNinspectedTVX{0};
  std::vector<std::string> swt_names;

  int mRunNumber;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before/", "after/"};

  void init(InitContext&)
  {
    mRunNumber = 0;
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    rctChecker.init(cfgRCTLabel.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);

    addhistograms();
  }

  ~eventQC() {}

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mTOIidx = zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfg_swt_name.value);
    for (auto& idx : mTOIidx) {
      LOGF(info, "Trigger of Interest : index = %d", idx);
    }
    mNinspectedTVX = zorro.getInspectedTVX()->GetBinContent(1);
    LOGF(info, "total inspected TVX events = %d in run number %d", mNinspectedTVX, bc.runNumber());
    fRegistry.fill(HIST("hNInspectedTVX"), bc.runNumber(), mNinspectedTVX);

    mRunNumber = bc.runNumber();
  }

  void addhistograms()
  {
    const AxisSpec axis_cent_ft0m{ConfCentBins, "centrality FT0M (%)"};
    const AxisSpec axis_cent_ft0a{ConfCentBins, "centrality FT0A (%)"};
    const AxisSpec axis_cent_ft0c{ConfCentBins, "centrality FT0C (%)"};

    const AxisSpec axis_mult_ft0a{ConfFT0AMultBins, "FT0A multiplicity"};
    const AxisSpec axis_mult_ft0c{ConfFT0CMultBins, "FT0C multiplicity"};
    const AxisSpec axis_mult_fv0a{ConfFV0AMultBins, "FV0A multiplicity"};
    const AxisSpec axis_mult_ncontrib{ConfTrackMultBins, "N_{track} to PV"};
    const AxisSpec axis_mult_ncontrib08{ConfTrackMultBins, "N_{track} to PV in |#eta| < 0.8"};
    const AxisSpec axis_mult_global_ncontrib08{ConfTrackMultBins, "N_{track}^{global} to PV in |#eta| < 0.8"};
    const AxisSpec axis_mult_globalTrack{ConfTrackMultBins, "N_{track}^{global} in |#eta| < 0.8"};

    if (doprocessEventQC_SWT) {
      fRegistry.add("BC/hNcoll", "Number of collisions per triggered BC;N_{collision} per triggered BC", kTH1F, {{11, -0.5, +10.5}}, false);
      fRegistry.add("BC/Collision/hMultNTracksPV", "hMultNTracksPV;N_{track} to PV", kTH1F, {{axis_mult_ncontrib}}, false);
      fRegistry.add("BC/Collision/hMultFT0AFT0C", "hMultFT0AFT0C;mult. FT0A;mult. FT0C", kTH2F, {{axis_mult_ft0a}, {axis_mult_ft0c}}, false);
      fRegistry.add("BC/Collision/hMultFT0AFV0A", "hMultFT0AFV0A;mult. FT0A;mult. FV0A", kTH2F, {{axis_mult_ft0a}, {axis_mult_fv0a}}, false);
      fRegistry.add("BC/Collision/hMultFT0CFV0A", "hMultFT0CFV0A;mult. FT0C;mult. FV0A", kTH2F, {{axis_mult_ft0c}, {axis_mult_fv0a}}, false);

      fRegistry.add("perBC/hDeltaTZ", "#DeltaZ_{vtx} vs. #DeltaT of collisions per BC;#DeltaZ_{vtx} (cm);#DeltaT (ns)", kTH2F, {{100, -5, +5}, {50, -25, +25}}, false);
      fRegistry.add("perBC/hCorrNcontrib", "hMultNTracksPV;", kTH2F, {{axis_mult_ncontrib}, {axis_mult_ncontrib}}, false);
      // fRegistry.addClone("perBC/", "beyondBC/");
    }

    // event info
    const int nbin_ev = 20;
    auto hCollisionCounter = fRegistry.add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{nbin_ev, 0.5, nbin_ev + 0.5}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "FT0AND");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No TF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No ITS ROF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "No Same Bunch Pileup");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "Is Vertex ITS-TPC");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "Is Vertex ITS-TPC-TRD");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "Is Vertex ITS-TPC-TOF");
    hCollisionCounter->GetXaxis()->SetBinLabel(10, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(11, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(12, "NoCollInTimeRangeStandard");
    hCollisionCounter->GetXaxis()->SetBinLabel(13, "NoCollInTimeRangeStrict");
    hCollisionCounter->GetXaxis()->SetBinLabel(14, "NoCollInRofStandard");
    hCollisionCounter->GetXaxis()->SetBinLabel(15, "NoCollInRofStrict");
    hCollisionCounter->GetXaxis()->SetBinLabel(16, "NoHighMultCollInPrevRof");
    hCollisionCounter->GetXaxis()->SetBinLabel(17, "GoodITSLayer3");
    hCollisionCounter->GetXaxis()->SetBinLabel(18, "GoodITSLayer0123");
    hCollisionCounter->GetXaxis()->SetBinLabel(19, "GoodITSLayersAll");
    hCollisionCounter->GetXaxis()->SetBinLabel(nbin_ev, "accepted");

    fRegistry.add("hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);

    if (cfgFillEvent) {
      fRegistry.add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
      fRegistry.add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{axis_mult_ncontrib}}, false);
      fRegistry.add("Event/before/hMultNTracksPV08", "hMultNTracksPV08; N_{track} to PV in |#eta| < 0.8", kTH1F, {{axis_mult_ncontrib08}}, false);
      fRegistry.add("Event/before/hMultFT0AFT0C", "hMultFT0AFT0C;mult. FT0A;mult. FT0C", kTH2F, {{axis_mult_ft0a}, {axis_mult_ft0c}}, false);
      fRegistry.add("Event/before/hMultFT0AFV0A", "hMultFT0AFV0A;mult. FT0A;mult. FV0A", kTH2F, {{axis_mult_ft0a}, {axis_mult_fv0a}}, false);
      fRegistry.add("Event/before/hMultFT0CFV0A", "hMultFT0CFV0A;mult. FT0C;mult. FV0A", kTH2F, {{axis_mult_ft0c}, {axis_mult_fv0a}}, false);
      fRegistry.add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{axis_cent_ft0a}}, false);
      fRegistry.add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{axis_cent_ft0c}}, false);
      fRegistry.add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{axis_cent_ft0m}}, false);
      fRegistry.add("Event/before/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV in |#eta| < 0.8", kTH2F, {{axis_cent_ft0c}, {axis_mult_ncontrib08}}, false);
      fRegistry.add("Event/before/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV in |#eta| < 0.8", kTH2F, {{axis_mult_ft0c}, {axis_mult_ncontrib08}}, false);
      fRegistry.add("Event/before/hMultFT0CvsTrackOccupancy", "hMultFT0CvsTrackOccupancy;mult. FT0C;N_{track} in time range", kTH2F, {{axis_mult_ft0c}, {200, 0, 20000}}, false);
      fRegistry.add("Event/before/hMultFV0AvsMultNTracksPV", "hMultFV0AvsMultNTracksPV;mult. FV0A;N_{track} to PV in |#eta| < 0.8", kTH2F, {{axis_mult_fv0a}, {axis_mult_ncontrib08}}, false);
      fRegistry.add("Event/before/hNTracksPVvsTrackOccupancy", "hNTracksPVvsTrackOccupancy;N_{track} to PV in |#eta| < 0.8;N_{track} in time range", kTH2F, {{axis_mult_ncontrib08}, {200, 0, 20000}}, false);
      fRegistry.add("Event/before/hNGlobalTracksvsTrackOccupancy", "hNGlobalTracksvsTrackOccupancy;N_{track}^{global} in |#eta| < 0.8;N_{track} in time range", kTH2F, {{axis_mult_globalTrack}, {200, 0, 20000}}, false);
      fRegistry.add("Event/before/hNGlobalTracksPVvsTrackOccupancy", "hNGlobalTracksPVvsTrackOccupancy;N_{track}^{global} to PV in |#eta| < 0.8;N_{track} in time range", kTH2F, {{axis_mult_global_ncontrib08}, {200, 0, 20000}}, false);
      fRegistry.add("Event/before/hCorrOccupancy", "occupancy correlation;FT0C occupancy;track-based occupancy", kTH2F, {{200, 0, 200000}, {200, 0, 20000}}, false);
    }
    fRegistry.addClone("Event/before/", "Event/after/");

    if (cfgFillEvent) {
      fRegistry.add("Event/after/hMultNGlobalTracks", "hMultNGlobalTracks; N_{track}^{global} in |#eta| < 0.8", kTH1F, {{axis_mult_globalTrack}}, false);
      fRegistry.add("Event/after/hCentFT0CvsMultNGlobalTracks", "hCentFT0CvsMultNGlobalTracks;centrality FT0C (%);N_{track}^{global} in |#eta| < 0.8", kTH2F, {{axis_cent_ft0c}, {axis_mult_globalTrack}}, false);
      fRegistry.add("Event/after/hMultFT0CvsMultNGlobalTracks", "hMultFT0CvsMultNGlobalTracks;mult. FT0C;N_{track}^{global} in |#eta| < 0.8", kTH2F, {{axis_mult_ft0c}, {axis_mult_globalTrack}}, false);
      fRegistry.add("Event/after/hMultNGlobalTracksPV", "hMultNGlobalTracksPV; N_{track}^{global} to PV in |#eta| < 0.8", kTH1F, {{axis_mult_global_ncontrib08}}, false);
      fRegistry.add("Event/after/hCentFT0CvsMultNGlobalTracksPV", "hCentFT0CvsMultNGlobalTracksPV;centrality FT0C (%);N_{track}^{global} to PV in |#eta| < 0.8", kTH2F, {{axis_cent_ft0c}, {axis_mult_global_ncontrib08}}, false);
      fRegistry.add("Event/after/hMultFT0CvsMultNGlobalTracksPV", "hMultFT0CvsMultNGlobalTracksPV;mult. FT0C;N_{track}^{global} to PV in |#eta| < 0.8", kTH2F, {{axis_mult_ft0c}, {axis_mult_global_ncontrib08}}, false);
    }

    std::vector<double> tmp_ptbins;
    for (int i = 0; i < 100; i++) {
      tmp_ptbins.emplace_back(0.01 * i); // every 0.01 GeV/c from 0 to 1 GeV/c
    }
    for (int i = 0; i < 91; i++) {
      tmp_ptbins.emplace_back(0.1 * i + 1.f); // every 0.1 GeV/c from 1 to 10 GeV/c
    }
    const AxisSpec axis_pt_tmp{tmp_ptbins, "p_{T} (GeV/c)"};

    const AxisSpec axis_pt{ConfPtBins, "p_{T} (GeV/c)"};
    const AxisSpec axis_eta{cfgNbinsEta, -2.0, +2.0, "#eta"};
    const AxisSpec axis_phi{cfgNbinsPhi, 0.0, 2 * M_PI, "#varphi (rad.)"};
    const AxisSpec axis_sign{3, -1.5, +1.5, "sign"};
    const AxisSpec axis_cent{20, 0, 100, "centrality FT0C (%)"};

    if (doprocessEventQC_Cent_Qvec) {
      for (int i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
        int nmod = cfgnMods->at(i);
        fRegistry.add(Form("Event/after/Qvector/hQ%dxFT0M_CentFT0C", nmod), Form("hQ%dxFT0M_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{FT0M}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dyFT0M_CentFT0C", nmod), Form("hQ%dyFT0M_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{FT0M}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dxFT0A_CentFT0C", nmod), Form("hQ%dxFT0A_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{FT0A}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dyFT0A_CentFT0C", nmod), Form("hQ%dyFT0A_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{FT0A}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dxFT0C_CentFT0C", nmod), Form("hQ%dxFT0C_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{FT0C}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dyFT0C_CentFT0C", nmod), Form("hQ%dyFT0C_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{FT0C}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dxFV0A_CentFT0C", nmod), Form("hQ%dxFV0A_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{FV0A}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dyFV0A_CentFT0C", nmod), Form("hQ%dyFV0A_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{FV0A}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dxBPos_CentFT0C", nmod), Form("hQ%dxBPos_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{BPos}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dyBPos_CentFT0C", nmod), Form("hQ%dyBPos_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{BPos}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dxBNeg_CentFT0C", nmod), Form("hQ%dxBNeg_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{BNeg}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dyBNeg_CentFT0C", nmod), Form("hQ%dyBNeg_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{BNeg}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dxBTot_CentFT0C", nmod), Form("hQ%dxBTot_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{BTot}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
        fRegistry.add(Form("Event/after/Qvector/hQ%dyBTot_CentFT0C", nmod), Form("hQ%dyBTot_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{BTot}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);

        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0MQ%dBPos_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0M} #upoint Q_{%d}^{BPos};centrality FT0C (%%);Q_{%d}^{FT0M} #upoint Q_{%d}^{BPos}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0MQ%dBNeg_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0M} #upoint Q_{%d}^{BNeg};centrality FT0C (%%);Q_{%d}^{FT0M} #upoint Q_{%d}^{BNeg}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dBPosQ%dBNeg_CentFT0C", nmod, nmod), Form("Q_{%d}^{BPos} #upoint Q_{%d}^{BNeg};centrality FT0C (%%);Q_{%d}^{BPos} #upoint Q_{%d}^{BNeg}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0CQ%dBPos_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0C} #upoint Q_{%d}^{BPos};centrality FT0C (%%);Q_{%d}^{FT0C} #upoint Q_{%d}^{BPos}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0CQ%dBNeg_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0C} #upoint Q_{%d}^{BNeg};centrality FT0C (%%);Q_{%d}^{FT0C} #upoint Q_{%d}^{BNeg}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0CQ%dBTot_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0C} #upoint Q_{%d}^{BTot};centrality FT0C (%%);Q_{%d}^{FT0C} #upoint Q_{%d}^{BTot}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0AQ%dBPos_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0A} #upoint Q_{%d}^{BPos};centrality FT0C (%%);Q_{%d}^{FT0A} #upoint Q_{%d}^{BPos}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0AQ%dBNeg_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0A} #upoint Q_{%d}^{BNeg};centrality FT0C (%%);Q_{%d}^{FT0A} #upoint Q_{%d}^{BNeg}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0AQ%dBTot_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0A} #upoint Q_{%d}^{BTot};centrality FT0C (%%);Q_{%d}^{FT0A} #upoint Q_{%d}^{BTot}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFT0AQ%dFT0C_CentFT0C", nmod, nmod), Form("Q_{%d}^{FT0A} #upoint Q_{%d}^{FT0C};centrality FT0C (%%);Q_{%d}^{FT0A} #upoint Q_{%d}^{FT0C}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFV0AQ%dBPos_CentFT0C", nmod, nmod), Form("Q_{%d}^{FV0A} #upoint Q_{%d}^{BPos};centrality FT0C (%%);Q_{%d}^{FV0A} #upoint Q_{%d}^{BPos}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFV0AQ%dBNeg_CentFT0C", nmod, nmod), Form("Q_{%d}^{FV0A} #upoint Q_{%d}^{BNeg};centrality FT0C (%%);Q_{%d}^{FV0A} #upoint Q_{%d}^{BNeg}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFV0AQ%dBTot_CentFT0C", nmod, nmod), Form("Q_{%d}^{FV0A} #upoint Q_{%d}^{BTot};centrality FT0C (%%);Q_{%d}^{FV0A} #upoint Q_{%d}^{BTot}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfQ%dFV0AQ%dFT0C_CentFT0C", nmod, nmod), Form("Q_{%d}^{FV0A} #upoint Q_{%d}^{FT0C};centrality FT0C (%%);Q_{%d}^{FV0A} #upoint Q_{%d}^{FT0C}", nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);

        fRegistry.add(Form("Event/after/Qvector/hEP%dFT0M_CentFT0C", nmod), Form("event plane FT0M;centrality FT0C (%%);#Psi_{%d}^{FT0M} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
        fRegistry.add(Form("Event/after/Qvector/hEP%dFT0A_CentFT0C", nmod), Form("event plane FT0A;centrality FT0C (%%);#Psi_{%d}^{FT0A} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
        fRegistry.add(Form("Event/after/Qvector/hEP%dFT0C_CentFT0C", nmod), Form("event plane FT0C;centrality FT0C (%%);#Psi_{%d}^{FT0C} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
        fRegistry.add(Form("Event/after/Qvector/hEP%dFV0A_CentFT0C", nmod), Form("event plane FV0A;centrality FT0C (%%);#Psi_{%d}^{FV0A} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
        fRegistry.add(Form("Event/after/Qvector/hEP%dBPos_CentFT0C", nmod), Form("event plane BPos;centrality FT0C (%%);#Psi_{%d}^{BPos} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
        fRegistry.add(Form("Event/after/Qvector/hEP%dBNeg_CentFT0C", nmod), Form("event plane BNeg;centrality FT0C (%%);#Psi_{%d}^{BNeg} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
        fRegistry.add(Form("Event/after/Qvector/hEP%dBTot_CentFT0C", nmod), Form("event plane BTot;centrality FT0C (%%);#Psi_{%d}^{BTot} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);

        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0MEP%dBPos_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0M} - #Psi_{%d}^{BPos}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0M} - #Psi_{%d}^{BPos}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0MEP%dBNeg_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0M} - #Psi_{%d}^{BNeg}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0M} - #Psi_{%d}^{BNeg}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dBPosEP%dBNeg_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{BPos} - #Psi_{%d}^{BNeg}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{BPos} - #Psi_{%d}^{BNeg}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0CEP%dBPos_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0C} - #Psi_{%d}^{BPos}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0C} - #Psi_{%d}^{BPos}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0CEP%dBNeg_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0C} - #Psi_{%d}^{BNeg}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0C} - #Psi_{%d}^{BNeg}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0CEP%dBTot_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0C} - #Psi_{%d}^{BTot}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0C} - #Psi_{%d}^{BTot}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0AEP%dBPos_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{BPos}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{BPos}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0AEP%dBNeg_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{BNeg}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{BNeg}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0AEP%dBTot_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{BTot}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{BTot}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFT0AEP%dFT0C_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{FT0C}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FT0A} - #Psi_{%d}^{FT0C}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFV0AEP%dBPos_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{BPos}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{BPos}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFV0AEP%dBNeg_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{BNeg}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{BNeg}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFV0AEP%dBTot_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{BTot}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{BTot}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
        fRegistry.add(Form("Event/after/Qvector/hPrfCosDiffEP%dFV0AEP%dFT0C_CentFT0C", nmod, nmod), Form("cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{FT0C}));centrality FT0C (%%);cos(%d(#Psi_{%d}^{FV0A} - #Psi_{%d}^{FT0C}))", nmod, nmod, nmod, nmod, nmod, nmod), kTProfile, {{100, 0, 100}}, false);
      }

      for (int i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
        int nmod = cfgnMods->at(i);
        const AxisSpec axis_sp{100, -5, +5, Form("#vec{u_{%d}} #upoint #vec{Q_{%d}}", nmod, nmod)};
        const AxisSpec axis_cos{200, -1, +1, Form("cos(%d(#varphi - #Psi_{%d}))", nmod, nmod)};
        fRegistry.add(Form("Track/hV%d", nmod), Form("charged particle v_{%d}", nmod), kTHnSparseD, {axis_cent, axis_pt, axis_sp, axis_cos}, false);
      }
    }

    if (cfgFillTrack) {
      fRegistry.add("Track/hs", "rec. single electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_sign}, false);
      fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -50, 50}}, false);
      fRegistry.add("Track/hRelSigma1Pt", "relative p_{T} resolution;p_{T} (GeV/c);#sigma_{1/p_{T}} #times p_{T}", kTH2F, {axis_pt_tmp, {100, 0, 0.1}}, false);
      fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {axis_pt_tmp, {500, 0., 500}}, false);
      fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {axis_pt_tmp, {500, 0., 500}}, false);
      fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hDeltaPin", "p_{in} vs. p_{pv};p_{pv} (GeV/c);(p_{in} - p_{pv})/p_{pv}", kTH2F, {{1000, 0, 10}, {200, -1, +1}}, false);
      fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("Track/hChi2TOF", "chi2 of TOF", kTH2F, {{1000, 0, 10}, {100, 0, 10}}, false);
    }

    if (cfgFillPID) {
      fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5.f, +5.f}}, false);
      fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5.f, +5.f}}, false);
      fRegistry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5.f, +5.f}}, false);
      fRegistry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5.f, +5.f}}, false);
      fRegistry.add("Track/hTOFbeta", "TOF #beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
      fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaKa", "TOF n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPr", "TOF n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

      fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda);", kTH2F, {{1000, 0.f, 10.f}, {150, 0, 15}}, false);
    }
  }

  template <typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    if (cfgFillTrack) {
      fRegistry.fill(HIST("Track/hs"), track.pt(), track.eta(), track.phi(), track.sign());
      fRegistry.fill(HIST("Track/hQoverPt"), track.signed1Pt());
      fRegistry.fill(HIST("Track/hRelSigma1Pt"), track.pt(), track.sigma1Pt() * track.pt());
      fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/hDCAxyzSigma"), track.dcaXY() / std::sqrt(track.cYY()), track.dcaZ() / std::sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/hDCAzRes_Pt"), track.pt(), std::sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
      fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/hDeltaPin"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
      fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/hChi2TOF"), track.p(), track.tofChi2());
    }
    if (cfgFillPID) {
      fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());

      fRegistry.fill(HIST("Track/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());

      int nsize = 0;
      for (int il = 0; il < 7; il++) {
        nsize += track.itsClsSizeInLayer(il);
      }
      fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.p(), static_cast<float>(nsize) / static_cast<float>(track.itsNCls()) * std::cos(std::atan(track.tgl())));
    }
  }

  template <int ev_id, typename TCollision>
  void fillEventInfo(TCollision const& collision)
  {
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
    }
    if (collision.sel8()) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 10.0);
    }
    if (std::fabs(collision.posZ()) < 10.0) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 11.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 12.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 13.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 14.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 15.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 16.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 17.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 18.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 19.0);
    }
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.numContrib());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV08"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0AFT0C"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0AFV0A"), collision.multFT0A(), collision.multFV0A());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CFV0A"), collision.multFT0C(), collision.multFV0A());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFV0AvsMultNTracksPV"), collision.multFV0A(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsTrackOccupancy"), collision.multFT0C(), collision.trackOccupancyInTimeRange());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hNTracksPVvsTrackOccupancy"), collision.multNTracksPV(), collision.trackOccupancyInTimeRange());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCorrOccupancy"), collision.ft0cOccupancyInTimeRange(), collision.trackOccupancyInTimeRange());

    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV());

    if constexpr (ev_id == 1 && std::is_same_v<std::decay_t<TCollision>, FilteredMyCollision_Qvec>) {
      if (std::find(cfgnMods->begin(), cfgnMods->end(), 2) != cfgnMods->end()) {
        fillQvectorInfo<ev_id, 2>(collision);
      }
      if (std::find(cfgnMods->begin(), cfgnMods->end(), 3) != cfgnMods->end()) {
        fillQvectorInfo<ev_id, 3>(collision);
      }
      if (std::find(cfgnMods->begin(), cfgnMods->end(), 4) != cfgnMods->end()) {
        fillQvectorInfo<ev_id, 4>(collision);
      }
    }
  }

  template <int ev_id, int nmod, typename TCollision>
  void fillQvectorInfo(TCollision const& collision)
  {
    int idx = std::distance(cfgnMods->begin(), std::find(cfgnMods->begin(), cfgnMods->end(), nmod));
    float qxft0m = collision.qvecFT0MReVec()[idx], qxft0a = collision.qvecFT0AReVec()[idx], qxft0c = collision.qvecFT0CReVec()[idx], qxfv0a = collision.qvecFV0AReVec()[idx], qxbpos = collision.qvecBPosReVec()[idx], qxbneg = collision.qvecBNegReVec()[idx], qxbtot = collision.qvecBTotReVec()[idx];
    float qyft0m = collision.qvecFT0MImVec()[idx], qyft0a = collision.qvecFT0AImVec()[idx], qyft0c = collision.qvecFT0CImVec()[idx], qyfv0a = collision.qvecFV0AImVec()[idx], qybpos = collision.qvecBPosImVec()[idx], qybneg = collision.qvecBNegImVec()[idx], qybtot = collision.qvecBTotImVec()[idx];
    std::array<float, 2> qft0m = {qxft0m, qyft0m};
    std::array<float, 2> qft0a = {qxft0a, qyft0a};
    std::array<float, 2> qft0c = {qxft0c, qyft0c};
    std::array<float, 2> qfv0a = {qxfv0a, qyfv0a};
    std::array<float, 2> qbpos = {qxbpos, qybpos};
    std::array<float, 2> qbneg = {qxbneg, qybneg};
    std::array<float, 2> qbtot = {qxbtot, qybtot};
    std::vector<std::array<float, 2>> qvectors = {qft0m, qft0a, qft0c, qfv0a, qbpos, qbneg, qbtot};

    if (!isGoodQvector(qvectors)) {
      return;
    }

    if constexpr (nmod == 2) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xFT0M_CentFT0C"), collision.centFT0C(), qxft0m);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yFT0M_CentFT0C"), collision.centFT0C(), qyft0m);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xFT0A_CentFT0C"), collision.centFT0C(), qxft0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yFT0A_CentFT0C"), collision.centFT0C(), qyft0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xFT0C_CentFT0C"), collision.centFT0C(), qxft0c);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yFT0C_CentFT0C"), collision.centFT0C(), qyft0c);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xFV0A_CentFT0C"), collision.centFT0C(), qxfv0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yFV0A_CentFT0C"), collision.centFT0C(), qyfv0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xBPos_CentFT0C"), collision.centFT0C(), qxbpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yBPos_CentFT0C"), collision.centFT0C(), qybpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xBNeg_CentFT0C"), collision.centFT0C(), qxbneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yBNeg_CentFT0C"), collision.centFT0C(), qybneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xBTot_CentFT0C"), collision.centFT0C(), qxbtot);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yBTot_CentFT0C"), collision.centFT0C(), qybtot);

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0MQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0MQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2BPosQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qbpos, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0CQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0CQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0CQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qft0c));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FV0AQ2FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qft0c));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FV0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FV0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FV0AQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qbtot));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2FT0M_CentFT0C"), collision.centFT0C(), getEP(qxft0m, qyft0m, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2FT0A_CentFT0C"), collision.centFT0C(), getEP(qxft0a, qyft0a, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2FT0C_CentFT0C"), collision.centFT0C(), getEP(qxft0c, qyft0c, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2FV0A_CentFT0C"), collision.centFT0C(), getEP(qxfv0a, qyfv0a, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2BPos_CentFT0C"), collision.centFT0C(), getEP(qxbpos, qybpos, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2BNeg_CentFT0C"), collision.centFT0C(), getEP(qxbneg, qybneg, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2BTot_CentFT0C"), collision.centFT0C(), getEP(qxbtot, qybtot, nmod));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0MEP2BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0MEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2BPosEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxbpos, qybpos, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0CEP2BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0CEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0CEP2BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2FT0C_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxft0c, qyft0c, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FV0AEP2FT0C_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxft0c, qyft0c, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FV0AEP2BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FV0AEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FV0AEP2BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxbtot, qybtot, nmod))));
    } else if constexpr (nmod == 3) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xFT0M_CentFT0C"), collision.centFT0C(), qxft0m);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yFT0M_CentFT0C"), collision.centFT0C(), qyft0m);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xFT0A_CentFT0C"), collision.centFT0C(), qxft0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yFT0A_CentFT0C"), collision.centFT0C(), qyft0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xFT0C_CentFT0C"), collision.centFT0C(), qxft0c);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yFT0C_CentFT0C"), collision.centFT0C(), qyft0c);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xFV0A_CentFT0C"), collision.centFT0C(), qxfv0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yFV0A_CentFT0C"), collision.centFT0C(), qyfv0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xBPos_CentFT0C"), collision.centFT0C(), qxbpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yBPos_CentFT0C"), collision.centFT0C(), qybpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xBNeg_CentFT0C"), collision.centFT0C(), qxbneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yBNeg_CentFT0C"), collision.centFT0C(), qybneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xBTot_CentFT0C"), collision.centFT0C(), qxbtot);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yBTot_CentFT0C"), collision.centFT0C(), qybtot);

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0MQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0MQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3BPosQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qbpos, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0CQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0CQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0CQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qft0c));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FV0AQ3FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qft0c));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FV0AQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FV0AQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FV0AQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qfv0a, qbtot));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3FT0M_CentFT0C"), collision.centFT0C(), getEP(qxft0m, qyft0m, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3FT0A_CentFT0C"), collision.centFT0C(), getEP(qxft0a, qyft0a, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3FT0C_CentFT0C"), collision.centFT0C(), getEP(qxft0c, qyft0c, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3FV0A_CentFT0C"), collision.centFT0C(), getEP(qxfv0a, qyfv0a, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3BPos_CentFT0C"), collision.centFT0C(), getEP(qxbpos, qybpos, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3BNeg_CentFT0C"), collision.centFT0C(), getEP(qxbneg, qybneg, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3BTot_CentFT0C"), collision.centFT0C(), getEP(qxbtot, qybtot, nmod));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0MEP3BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0MEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3BPosEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxbpos, qybpos, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0CEP3BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0CEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0CEP3BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3FT0C_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxft0c, qyft0c, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FV0AEP3FT0C_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxft0c, qyft0c, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FV0AEP3BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FV0AEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FV0AEP3BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxfv0a, qyfv0a, nmod) - getEP(qxbtot, qybtot, nmod))));
    }
  }

  template <int nmod, typename TCollision, typename TTrack>
  void fillVn(TCollision const& collision, TTrack const& track)
  {
    int idx = std::distance(cfgnMods->begin(), std::find(cfgnMods->begin(), cfgnMods->end(), nmod));
    float qxft0m = collision.qvecFT0MReVec()[idx], qxft0a = collision.qvecFT0AReVec()[idx], qxft0c = collision.qvecFT0CReVec()[idx], qxfv0a = collision.qvecFV0AReVec()[idx], qxbpos = collision.qvecBPosReVec()[idx], qxbneg = collision.qvecBNegReVec()[idx], qxbtot = collision.qvecBTotReVec()[idx];
    float qyft0m = collision.qvecFT0MImVec()[idx], qyft0a = collision.qvecFT0AImVec()[idx], qyft0c = collision.qvecFT0CImVec()[idx], qyfv0a = collision.qvecFV0AImVec()[idx], qybpos = collision.qvecBPosImVec()[idx], qybneg = collision.qvecBNegImVec()[idx], qybtot = collision.qvecBTotImVec()[idx];
    std::array<float, 2> qft0m = {qxft0m, qyft0m};
    std::array<float, 2> qft0a = {qxft0a, qyft0a};
    std::array<float, 2> qft0c = {qxft0c, qyft0c};
    std::array<float, 2> qfv0a = {qxfv0a, qyfv0a};
    std::array<float, 2> qbpos = {qxbpos, qybpos};
    std::array<float, 2> qbneg = {qxbneg, qybneg};
    std::array<float, 2> qbtot = {qxbtot, qybtot};
    std::vector<std::array<float, 2>> qvectors = {qft0m, qft0a, qft0c, qfv0a, qbpos, qbneg, qbtot};

    if (!isGoodQvector(qvectors)) {
      return;
    }

    const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};

    float sp = RecoDecay::dotProd(std::array<float, 2>{static_cast<float>(std::cos(nmod * track.phi())), static_cast<float>(std::sin(nmod * track.phi()))}, qvectors[cfgQvecEstimator]);
    float cos = std::cos(nmod * (track.phi() - getEP(qvectors[cfgQvecEstimator][0], qvectors[cfgQvecEstimator][1], nmod)));
    if constexpr (nmod == 2) {
      fRegistry.fill(HIST("Track/hV2"), centralities[cfgCentEstimator], track.pt(), sp, cos);
    } else if constexpr (nmod == 3) {
      fRegistry.fill(HIST("Track/hV3"), centralities[cfgCentEstimator], track.pt(), sp, cos);
    }
  }

  float getEP(float qx, float qy, int nmod)
  {
    return std::atan2(qy, qx) / static_cast<float>(nmod);
  }

  template <typename TQvectors>
  bool isGoodQvector(TQvectors const& qvectors)
  {
    bool is_good = true;
    for (auto& qvec : qvectors) {
      if (std::fabs(qvec[0]) > 20.f || std::fabs(qvec[1]) > 20.f) {
        is_good = false;
        break;
      }
    }
    return is_good;
  }

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.pt() < trackcuts.cfg_min_pt_track || trackcuts.cfg_max_pt_track < track.pt()) {
      return false;
    }

    if (track.eta() < trackcuts.cfg_min_eta_track || trackcuts.cfg_max_eta_track < track.eta()) {
      return false;
    }

    if (std::fabs(track.dcaXY()) > trackcuts.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(track.dcaZ()) > trackcuts.cfg_max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() > trackcuts.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < trackcuts.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < trackcuts.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > trackcuts.cfg_max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < trackcuts.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < trackcuts.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < trackcuts.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > trackcuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    if (track.tpcNSigmaEl() < trackcuts.cfg_min_TPCNsigmaEl || trackcuts.cfg_max_TPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }

    if (trackcuts.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < trackcuts.cfg_max_TPCNsigmaPi) {
      return false;
    }

    if (trackcuts.cfg_requireTOF && !(track.hasTOF() && track.tofChi2() < trackcuts.cfg_max_chi2tof)) {
      return false;
    }

    if (track.hasTOF() && ((track.tofNSigmaEl() < trackcuts.cfg_min_TOFNsigmaEl || trackcuts.cfg_max_TOFNsigmaEl < track.tofNSigmaEl()) || trackcuts.cfg_max_chi2tof < track.tofChi2())) {
      return false;
    }

    return true;
  }

  template <typename TCollision>
  bool isSelectedCollision(TCollision const& collision)
  {
    if (eventcuts.cfgRequireSel8 && !collision.sel8()) {
      return false;
    }

    if (collision.posZ() < eventcuts.cfgZvtxMin || eventcuts.cfgZvtxMax < collision.posZ()) {
      return false;
    }

    if (eventcuts.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (eventcuts.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (eventcuts.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (eventcuts.cfgRequireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }

    if (eventcuts.cfgRequireVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }

    if (eventcuts.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (eventcuts.cfgRequireNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }

    if (eventcuts.cfgRequireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }

    if (eventcuts.cfgRequirekNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }

    if (eventcuts.cfgRequirekNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }

    if (eventcuts.cfgRequirekNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodITSLayer3 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodITSLayer0123 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }

    if (!(eventcuts.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax)) {
      return false;
    }

    if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      return false;
    }

    return true;
  }

  Filter collisionFilter_evsel = ifnode(eventcuts.cfgRequireSel8.node(), o2::aod::evsel::sel8 == true, true);
  Filter collisionFilter_zvtx = eventcuts.cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < eventcuts.cfgZvtxMax;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_numContrib = cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < cfgNumContribMax;
  Filter collisionFilter_track_occupancy = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_ft0c_occupancy = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;
  using FilteredMyCollisions_Qvec = soa::Filtered<MyCollisions_Qvec>;

  using FilteredMyCollision = FilteredMyCollisions::iterator;
  using FilteredMyCollision_Qvec = FilteredMyCollisions_Qvec::iterator;

  Filter trackFilter = (trackcuts.cfg_min_pt_track < o2::aod::track::pt && o2::aod::track::pt < trackcuts.cfg_max_pt_track) && (trackcuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < trackcuts.cfg_max_eta_track) && nabs(o2::aod::track::dcaXY) < trackcuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < trackcuts.cfg_max_dcaz && o2::aod::track::tpcChi2NCl < trackcuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < trackcuts.cfg_max_chi2its && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  using FilteredMyTracks = soa::Filtered<MyTracks>;

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  // Preslice<aod::Collisions> perBC = o2::aod::collision::bcId;
  PresliceUnsorted<MyCollisions> perFoundBC = aod::evsel::foundBCId;

  template <bool isTriggerAnalysis, typename TBCs, typename TCollisions, typename TTracks>
  void runQC(TBCs const& bcs, TCollisions const& collisions, TTracks const& tracks)
  {
    if constexpr (isTriggerAnalysis) {
      // std::vector<int> selectedCollisionIds;
      // selectedCollisionIds.reserve(collisions.size());

      for (const auto& bc : bcs) {
        initCCDB(bc);
        if (!zorro.isSelected(bc.globalBC())) { // triggered BC
          continue;
        }

        // const auto& collisions_per_bc = collisions.sliceBy(perBC, bc.globalIndex());
        const auto& collisions_per_bc = collisions.sliceBy(perFoundBC, bc.globalIndex());
        fRegistry.fill(HIST("BC/hNcoll"), collisions_per_bc.size());
        for (const auto& collision : collisions_per_bc) {
          if (!isSelectedCollision(collision)) {
            continue;
          }

          fRegistry.fill(HIST("BC/Collision/hMultNTracksPV"), collision.numContrib());
          fRegistry.fill(HIST("BC/Collision/hMultFT0AFT0C"), collision.multFT0A(), collision.multFT0C());
          fRegistry.fill(HIST("BC/Collision/hMultFT0AFV0A"), collision.multFT0A(), collision.multFV0A());
          fRegistry.fill(HIST("BC/Collision/hMultFT0CFV0A"), collision.multFT0C(), collision.multFV0A());
          // selectedCollisionIds.emplace_back(collision.globalIndex());
        }

        for (const auto& [col1, col2] : combinations(CombinationsStrictlyUpperIndexPolicy(collisions_per_bc, collisions_per_bc))) {
          if (!isSelectedCollision(col1) || !isSelectedCollision(col2)) {
            continue;
          }
          fRegistry.fill(HIST("perBC/hDeltaTZ"), col1.posZ() - col2.posZ(), col1.collisionTime() - col2.collisionTime());
          fRegistry.fill(HIST("perBC/hCorrNcontrib"), col1.numContrib(), col2.numContrib());
        } // end of pairing
      } // end of bc loop

      // for (const auto& collisionId1 : selectedCollisionIds) {
      //   const auto& col1 = collisions.rawIteratorAt(collisionId1);
      //   for (const auto& col2 : collisions) {
      //     if (!isSelectedCollision(col2)) {
      //       continue;
      //     }

      //     const auto& bc1 = col1.template bc_as<TBCs>(); // don't use foundBC for CEFP.
      //     const auto& bc2 = col2.template bc_as<TBCs>(); // don't use foundBC for CEFP.
      //     if (bc1.globalBC() == bc2.globalBC()) {
      //       continue;
      //     }
      //     fRegistry.fill(HIST("beyondBC/hDeltaTZ"), col1.posZ() - col2.posZ(), col1.collisionTime() - col2.collisionTime());
      //     fRegistry.fill(HIST("beyondBC/hCorrNcontrib"), col1.numContrib(), col2.numContrib());
      //   } // end of all collision loop
      // } // end of selected collision loop

      // selectedCollisionIds.clear();
      // selectedCollisionIds.shrink_to_fit();
    } // end of trigger QC

    for (const auto& collision : collisions) {
      if constexpr (isTriggerAnalysis) {
        const auto& bc = collision.template bc_as<TBCs>(); // don't use foundBC for CEFP.
        initCCDB(bc);
        if (!zorro.isSelected(bc.globalBC())) { // triggered event
          continue;
        }
      }

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      if (cfgFillEvent) {
        fillEventInfo<0>(collision);
      }
      if (!isSelectedCollision(collision)) {
        continue;
      }
      if (cfgFillEvent) {
        fillEventInfo<1>(collision);
      }
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 20); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 20);  // accepted

      int nGlobalTracks = 0, nGlobalTracksPV = 0;
      auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
      for (const auto& track : tracks_per_coll) {
        if (!isSelectedTrack(track)) {
          continue;
        }
        if (isElectron(track)) {
          fillTrackInfo(track);
        }
        if (std::fabs(track.eta()) < 0.8) {
          nGlobalTracks++;

          if constexpr (std::is_same_v<std::decay_t<TCollisions>, FilteredMyCollisions_Qvec>) {
            for (int i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
              if (cfgnMods->at(i) == 2) {
                fillVn<2>(collision, track);
              } else if (cfgnMods->at(i) == 3) {
                fillVn<3>(collision, track);
              }
            }
          }

          if (track.isPVContributor()) {
            nGlobalTracksPV++;
          }
        }
      } // end of track loop

      if (cfgFillEvent) {
        fRegistry.fill(HIST("Event/after/hMultNGlobalTracks"), nGlobalTracks);
        fRegistry.fill(HIST("Event/after/hMultNGlobalTracksPV"), nGlobalTracksPV);
        fRegistry.fill(HIST("Event/after/hMultFT0CvsMultNGlobalTracks"), collision.multFT0C(), nGlobalTracks);
        fRegistry.fill(HIST("Event/after/hMultFT0CvsMultNGlobalTracksPV"), collision.multFT0C(), nGlobalTracksPV);
        fRegistry.fill(HIST("Event/after/hNGlobalTracksvsTrackOccupancy"), nGlobalTracks, collision.trackOccupancyInTimeRange());
        fRegistry.fill(HIST("Event/after/hNGlobalTracksPVvsTrackOccupancy"), nGlobalTracksPV, collision.trackOccupancyInTimeRange());
        fRegistry.fill(HIST("Event/after/hCentFT0CvsMultNGlobalTracks"), collision.centFT0C(), nGlobalTracks);
        fRegistry.fill(HIST("Event/after/hCentFT0CvsMultNGlobalTracksPV"), collision.centFT0C(), nGlobalTracksPV);
      }

    } // end of collision loop
  } // end of process

  void processEventQC(MyBCs const& bcs, FilteredMyCollisions const& collisions, FilteredMyTracks const& tracks)
  {
    runQC<false>(bcs, collisions, tracks);
  }
  PROCESS_SWITCH(eventQC, processEventQC, "event QC", true);

  void processEventQC_SWT(MyBCs const& bcs, FilteredMyCollisions const& collisions, FilteredMyTracks const& tracks)
  {
    runQC<true>(bcs, collisions, tracks);
  }
  PROCESS_SWITCH(eventQC, processEventQC_SWT, "event QC", false);

  void processEventQC_Cent_Qvec(MyBCs const& bcs, FilteredMyCollisions_Qvec const& collisions, FilteredMyTracks const& tracks)
  {
    runQC<false>(bcs, collisions, tracks);
  }
  PROCESS_SWITCH(eventQC, processEventQC_Cent_Qvec, "event QC + q vector", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<eventQC>(cfgc, TaskName{"event-qc"})};
}
