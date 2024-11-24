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

#include <string>
#include <vector>
#include <algorithm>
#include <type_traits>

#include "TString.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "MathUtils/Utils.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
using MyCollisions_Qvec = soa::Join<MyCollisions, MyQvectors>;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;

struct eventQC {
  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<bool> cfgFillPID{"cfgFillPID", false, "fill PID histograms"};
  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2, 3}, "Modulation of interest. Please keep increasing order"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  ConfigurableAxis ConfPtBins{"ConfPtBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pT bins for output histograms"};
  Configurable<int> cfgNbinsEta{"cfgNbinsEta", 40, "number of eta bins for output histograms"};
  Configurable<int> cfgNbinsPhi{"cfgNbinsPhi", 360, "number of phi bins for output histograms"};

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
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
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster its"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
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

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before/", "after/"};

  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    addhistograms();
  }

  ~eventQC() {}

  void addhistograms()
  {
    // event info

    const int nbin_ev = 13;
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
    hCollisionCounter->GetXaxis()->SetBinLabel(13, "accepted");

    fRegistry.add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{200, 0, 200000}, {60, 0, 60000}}, false);
    fRegistry.add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", kTH2F, {{100, 0, 100}, {600, 0, 6000}}, false);
    fRegistry.add("Event/before/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", kTH2F, {{60, 0, 60000}, {600, 0, 6000}}, false);
    fRegistry.add("Event/before/hMultFT0CvsOccupancy", "hMultFT0CvsOccupancy;mult. FT0C;N_{track} in time range", kTH2F, {{60, 0, 60000}, {200, 0, 20000}}, false);
    fRegistry.add("Event/before/hNTracksPVvsOccupancy", "hNTracksPVvsOccupancy;N_{track} to PV;N_{track} in time range", kTH2F, {{600, 0, 6000}, {200, 0, 20000}}, false);
    fRegistry.add("Event/before/hNGlobalTracksvsOccupancy", "hNGlobalTracksvsOccupancy;N_{track}^{global};N_{track} in time range", kTH2F, {{600, 0, 6000}, {200, 0, 20000}}, false);
    fRegistry.add("Event/before/hNGlobalTracksPVvsOccupancy", "hNGlobalTracksPVvsOccupancy;N_{track}^{global} to PV;N_{track} in time range", kTH2F, {{600, 0, 6000}, {200, 0, 20000}}, false);
    fRegistry.add("Event/before/hCorrOccupancy", "occupancy correlation;FT0C occupancy;track-based occupancy", kTH2F, {{200, 0, 200000}, {200, 0, 20000}}, false);
    fRegistry.addClone("Event/before/", "Event/after/");

    fRegistry.add("Event/after/hMultNGlobalTracks", "hMultNGlobalTracks; N_{track}^{global}", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/after/hCentFT0CvsMultNGlobalTracks", "hCentFT0CvsMultNGlobalTracks;centrality FT0C (%);N_{track}^{global}", kTH2F, {{100, 0, 100}, {600, 0, 6000}}, false);
    fRegistry.add("Event/after/hMultFT0CvsMultNGlobalTracks", "hMultFT0CvsMultNGlobalTracks;mult. FT0C;N_{track}^{global}", kTH2F, {{60, 0, 60000}, {600, 0, 6000}}, false);
    fRegistry.add("Event/after/hMultNGlobalTracksPV", "hMultNGlobalTracksPV; N_{track}^{global} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/after/hCentFT0CvsMultNGlobalTracksPV", "hCentFT0CvsMultNGlobalTracksPV;centrality FT0C (%);N_{track}^{global} to PV", kTH2F, {{100, 0, 100}, {600, 0, 6000}}, false);
    fRegistry.add("Event/after/hMultFT0CvsMultNGlobalTracksPV", "hMultFT0CvsMultNGlobalTracksPV;mult. FT0C;N_{track}^{global} to PV", kTH2F, {{60, 0, 60000}, {600, 0, 6000}}, false);

    for (int i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
      int nmod = cfgnMods->at(i);
      fRegistry.add(Form("Event/after/Qvector/hQ%dxFT0M_CentFT0C", nmod), Form("hQ%dxFT0M_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{FT0M}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
      fRegistry.add(Form("Event/after/Qvector/hQ%dyFT0M_CentFT0C", nmod), Form("hQ%dyFT0M_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{FT0M}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
      fRegistry.add(Form("Event/after/Qvector/hQ%dxFT0A_CentFT0C", nmod), Form("hQ%dxFT0A_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{FT0A}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
      fRegistry.add(Form("Event/after/Qvector/hQ%dyFT0A_CentFT0C", nmod), Form("hQ%dyFT0A_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{FT0A}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
      fRegistry.add(Form("Event/after/Qvector/hQ%dxFT0C_CentFT0C", nmod), Form("hQ%dxFT0C_CentFT0C;centrality FT0C (%%);Q_{%d,x}^{FT0C}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
      fRegistry.add(Form("Event/after/Qvector/hQ%dyFT0C_CentFT0C", nmod), Form("hQ%dyFT0C_CentFT0C;centrality FT0C (%%);Q_{%d,y}^{FT0C}", nmod, nmod), kTH2F, {{100, 0, 100}, {200, -10, +10}}, false);
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

      fRegistry.add(Form("Event/after/Qvector/hEP%dFT0M_CentFT0C", nmod), Form("event plane FT0M;centrality FT0C (%%);#Psi_{%d}^{FT0M} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
      fRegistry.add(Form("Event/after/Qvector/hEP%dFT0A_CentFT0C", nmod), Form("event plane FT0A;centrality FT0C (%%);#Psi_{%d}^{FT0A} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
      fRegistry.add(Form("Event/after/Qvector/hEP%dFT0C_CentFT0C", nmod), Form("event plane FT0C;centrality FT0C (%%);#Psi_{%d}^{FT0C} (rad.)", nmod), kTH2F, {{100, 0, 100}, {36, -M_PI_2, +M_PI_2}}, false);
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
    const AxisSpec axis_eta{cfgNbinsEta, -1.0, +1.0, "#eta"};
    const AxisSpec axis_phi{cfgNbinsPhi, 0.0, 2 * M_PI, "#varphi (rad.)"};
    const AxisSpec axis_sign{3, -1.5, +1.5, "sign"};
    const AxisSpec axis_cent{20, 0, 100, "centrality FT0C (%)"};

    for (int i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
      int nmod = cfgnMods->at(i);
      const AxisSpec axis_sp{100, -5, +5, Form("#vec{u_{%d}} #upoint #vec{Q_{%d}}", nmod, nmod)};
      const AxisSpec axis_cos{200, -1, +1, Form("cos(%d(#varphi - #Psi_{%d}))", nmod, nmod)};
      fRegistry.add(Form("Track/hV%d", nmod), Form("charged particle v_{%d}", nmod), kTHnSparseD, {axis_cent, axis_pt, axis_sp, axis_cos}, false);
    }
    fRegistry.add("Track/hs", "rec. single electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_sign}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
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
    fRegistry.add("Track/hChi2TOF", "chi2 of TOF", kTH1F, {{100, 0, 10}}, false);

    if (cfgFillPID) {
      fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{pv} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/hTOFbeta", "TOF #beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
      fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda);", kTH2F, {{1000, 0.f, 10.f}, {150, 0, 15}}, false);
      fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaMu", "TPC n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaMu", "TOF n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaKa", "TOF n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPr", "TOF n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    }
  }

  template <typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    fRegistry.fill(HIST("Track/hs"), track.pt(), track.eta(), track.phi(), track.sign());
    fRegistry.fill(HIST("Track/hQoverPt"), track.signed1Pt());
    fRegistry.fill(HIST("Track/hRelSigma1Pt"), track.pt(), track.sigma1Pt() * track.pt());
    fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
    fRegistry.fill(HIST("Track/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
    fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
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
    fRegistry.fill(HIST("Track/hChi2TOF"), track.tofChi2());

    if (cfgFillPID) {
      int nsize = 0;
      for (int il = 0; il < 7; il++) {
        nsize += track.itsClsSizeInLayer(il);
      }
      fRegistry.fill(HIST("Track/hTPCdEdx"), track.p(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.p(), static_cast<float>(nsize) / static_cast<float>(track.itsNCls()) * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.p(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTPCNsigmaMu"), track.p(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.p(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/hTPCNsigmaKa"), track.p(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/hTPCNsigmaPr"), track.p(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
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
    if (fabs(collision.posZ()) < 10.0) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 11.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 12.0);
    }
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0CvsOccupancy"), collision.multFT0C(), collision.trackOccupancyInTimeRange());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hNTracksPVvsOccupancy"), collision.multNTracksPV(), collision.trackOccupancyInTimeRange());
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
    float qxft0m = collision.qvecFT0MReVec()[idx], qxft0a = collision.qvecFT0AReVec()[idx], qxft0c = collision.qvecFT0CReVec()[idx], qxbpos = collision.qvecBPosReVec()[idx], qxbneg = collision.qvecBNegReVec()[idx], qxbtot = collision.qvecBTotReVec()[idx];
    float qyft0m = collision.qvecFT0MImVec()[idx], qyft0a = collision.qvecFT0AImVec()[idx], qyft0c = collision.qvecFT0CImVec()[idx], qybpos = collision.qvecBPosImVec()[idx], qybneg = collision.qvecBNegImVec()[idx], qybtot = collision.qvecBTotImVec()[idx];
    std::array<float, 2> qft0m = {qxft0m, qyft0m};
    std::array<float, 2> qft0a = {qxft0a, qyft0a};
    std::array<float, 2> qft0c = {qxft0c, qyft0c};
    std::array<float, 2> qbpos = {qxbpos, qybpos};
    std::array<float, 2> qbneg = {qxbneg, qybneg};
    std::array<float, 2> qbtot = {qxbtot, qybtot};
    std::vector<std::array<float, 2>> qvectors = {qft0m, qft0a, qft0c, qbpos, qbneg, qbtot};

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
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xBPos_CentFT0C"), collision.centFT0C(), qxbpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yBPos_CentFT0C"), collision.centFT0C(), qybpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xBNeg_CentFT0C"), collision.centFT0C(), qxbneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yBNeg_CentFT0C"), collision.centFT0C(), qybneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2xBTot_CentFT0C"), collision.centFT0C(), qxbtot);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ2yBTot_CentFT0C"), collision.centFT0C(), qybtot);

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0MQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0MQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2BPosQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qbpos, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0CQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0CQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0CQ2BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ2FT0AQ2FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qft0c));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2FT0M_CentFT0C"), collision.centFT0C(), getEP(qxft0m, qyft0m, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2FT0A_CentFT0C"), collision.centFT0C(), getEP(qxft0a, qyft0a, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2FT0C_CentFT0C"), collision.centFT0C(), getEP(qxft0c, qyft0c, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2BPos_CentFT0C"), collision.centFT0C(), getEP(qxbpos, qybpos, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2BNeg_CentFT0C"), collision.centFT0C(), getEP(qxbneg, qybneg, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP2BTot_CentFT0C"), collision.centFT0C(), getEP(qxbtot, qybtot, nmod));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0MEP2BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0MEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2BPosEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxbpos, qybpos, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0CEP2BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0CEP2BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0CEP2BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP2FT0AEP2FT0C_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxft0c, qyft0c, nmod))));
    } else if constexpr (nmod == 3) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xFT0M_CentFT0C"), collision.centFT0C(), qxft0m);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yFT0M_CentFT0C"), collision.centFT0C(), qyft0m);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xFT0A_CentFT0C"), collision.centFT0C(), qxft0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yFT0A_CentFT0C"), collision.centFT0C(), qyft0a);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xFT0C_CentFT0C"), collision.centFT0C(), qxft0c);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yFT0C_CentFT0C"), collision.centFT0C(), qyft0c);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xBPos_CentFT0C"), collision.centFT0C(), qxbpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yBPos_CentFT0C"), collision.centFT0C(), qybpos);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xBNeg_CentFT0C"), collision.centFT0C(), qxbneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yBNeg_CentFT0C"), collision.centFT0C(), qybneg);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3xBTot_CentFT0C"), collision.centFT0C(), qxbtot);
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hQ3yBTot_CentFT0C"), collision.centFT0C(), qybtot);

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0MQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0MQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0m, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3BPosQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qbpos, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0CQ3BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbpos));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0CQ3BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbneg));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0CQ3BTot_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0c, qbtot));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfQ3FT0AQ3FT0C_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(qft0a, qft0c));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3FT0M_CentFT0C"), collision.centFT0C(), getEP(qxft0m, qyft0m, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3FT0A_CentFT0C"), collision.centFT0C(), getEP(qxft0a, qyft0a, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3FT0C_CentFT0C"), collision.centFT0C(), getEP(qxft0c, qyft0c, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3BPos_CentFT0C"), collision.centFT0C(), getEP(qxbpos, qybpos, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3BNeg_CentFT0C"), collision.centFT0C(), getEP(qxbneg, qybneg, nmod));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hEP3BTot_CentFT0C"), collision.centFT0C(), getEP(qxbtot, qybtot, nmod));

      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0MEP3BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0MEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0m, qyft0m, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3BPosEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxbpos, qybpos, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0CEP3BPos_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbpos, qybpos, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0CEP3BNeg_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbneg, qybneg, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0CEP3BTot_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0c, qyft0c, nmod) - getEP(qxbtot, qybtot, nmod))));
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("Qvector/hPrfCosDiffEP3FT0AEP3FT0C_CentFT0C"), collision.centFT0C(), std::cos(nmod * (getEP(qxft0a, qyft0a, nmod) - getEP(qxft0c, qyft0c, nmod))));
    }
  }

  template <int nmod, typename TCollision, typename TTrack>
  void fillVn(TCollision const& collision, TTrack const& track)
  {
    int idx = std::distance(cfgnMods->begin(), std::find(cfgnMods->begin(), cfgnMods->end(), nmod));
    float qxft0m = collision.qvecFT0MReVec()[idx], qxft0a = collision.qvecFT0AReVec()[idx], qxft0c = collision.qvecFT0CReVec()[idx], qxbpos = collision.qvecBPosReVec()[idx], qxbneg = collision.qvecBNegReVec()[idx], qxbtot = collision.qvecBTotReVec()[idx];
    float qyft0m = collision.qvecFT0MImVec()[idx], qyft0a = collision.qvecFT0AImVec()[idx], qyft0c = collision.qvecFT0CImVec()[idx], qybpos = collision.qvecBPosImVec()[idx], qybneg = collision.qvecBNegImVec()[idx], qybtot = collision.qvecBTotImVec()[idx];
    std::array<float, 2> qft0m = {qxft0m, qyft0m};
    std::array<float, 2> qft0a = {qxft0a, qyft0a};
    std::array<float, 2> qft0c = {qxft0c, qyft0c};
    std::array<float, 2> qbpos = {qxbpos, qybpos};
    std::array<float, 2> qbneg = {qxbneg, qybneg};
    std::array<float, 2> qbtot = {qxbtot, qybtot};
    std::vector<std::array<float, 2>> qvectors = {qft0m, qft0a, qft0c, qbpos, qbneg, qbtot};

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
      if (fabs(qvec[0]) > 20.f || fabs(qvec[1]) > 20.f) {
        is_good = false;
        break;
      }
    }
    return is_good;
  }

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track)
  {
    if (track.itsNCls() < trackcuts.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < trackcuts.cfg_min_ncluster_itsib) {
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

    if (trackcuts.cfg_requireTOF && !track.hasTOF()) {
      return false;
    }

    if (track.hasTOF() && (track.tofNSigmaEl() < trackcuts.cfg_min_TOFNsigmaEl || trackcuts.cfg_max_TOFNsigmaEl < track.tofNSigmaEl())) {
      return false;
    }

    return true;
  }

  template <typename TCollision>
  bool isSelectedEvent(TCollision const& collision)
  {
    if (!collision.sel8()) {
      return false;
    }

    if (fabs(collision.posZ()) > eventcuts.cfgZvtxMax) {
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

    if (!(eventcuts.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax)) {
      return false;
    }

    if (!(eventcuts.cfgFT0COccupancyMin < collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      return false;
    }

    return true;
  }

  Filter collisionFilter_evsel = o2::aod::evsel::sel8 == true && nabs(o2::aod::collision::posZ) < eventcuts.cfgZvtxMax;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_track_occupancy = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_ft0c_occupancy = eventcuts.cfgFT0COccupancyMin < o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;
  using FilteredMyCollisions_Qvec = soa::Filtered<MyCollisions_Qvec>;

  using FilteredMyCollision = FilteredMyCollisions::iterator;
  using FilteredMyCollision_Qvec = FilteredMyCollisions_Qvec::iterator;

  Filter trackFilter = (trackcuts.cfg_min_pt_track < o2::aod::track::pt && o2::aod::track::pt < trackcuts.cfg_max_pt_track) && (trackcuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < trackcuts.cfg_max_eta_track) && nabs(o2::aod::track::dcaXY) < trackcuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < trackcuts.cfg_max_dcaz && o2::aod::track::tpcChi2NCl < trackcuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < trackcuts.cfg_max_chi2its && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  using FilteredMyTracks = soa::Filtered<MyTracks>;

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;

  template <typename TCollisions>
  void processQC(TCollisions const& collisions, FilteredMyTracks const& tracks)
  {
    for (auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fillEventInfo<0>(collision);
      if (!isSelectedEvent(collision)) {
        continue;
      }
      fillEventInfo<1>(collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 13); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 13);  // accepted

      int nGlobalTracks = 0, nGlobalTracksPV = 0;
      auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
      for (auto& track : tracks_per_coll) {
        if (!isSelectedTrack(track)) {
          continue;
        }
        if (isElectron(track)) {
          fillTrackInfo(track);
        }
        if (fabs(track.eta()) < 0.8) {
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
      }

      fRegistry.fill(HIST("Event/after/hMultNGlobalTracks"), nGlobalTracks);
      fRegistry.fill(HIST("Event/after/hMultNGlobalTracksPV"), nGlobalTracksPV);
      fRegistry.fill(HIST("Event/after/hMultFT0CvsMultNGlobalTracks"), collision.multFT0C(), nGlobalTracks);
      fRegistry.fill(HIST("Event/after/hMultFT0CvsMultNGlobalTracksPV"), collision.multFT0C(), nGlobalTracksPV);
      fRegistry.fill(HIST("Event/after/hNGlobalTracksvsOccupancy"), nGlobalTracks, collision.trackOccupancyInTimeRange());
      fRegistry.fill(HIST("Event/after/hNGlobalTracksPVvsOccupancy"), nGlobalTracksPV, collision.trackOccupancyInTimeRange());
      if constexpr (std::is_same_v<std::decay_t<TCollisions>, FilteredMyCollisions_Qvec>) {
        fRegistry.fill(HIST("Event/after/hCentFT0CvsMultNGlobalTracks"), collision.centFT0C(), nGlobalTracks);
        fRegistry.fill(HIST("Event/after/hCentFT0CvsMultNGlobalTracksPV"), collision.centFT0C(), nGlobalTracksPV);
      }
    } // end of collision loop
  } // end of process

  PROCESS_SWITCH_FULL(eventQC, processQC<FilteredMyCollisions>, processEventQC, "event QC", true);
  PROCESS_SWITCH_FULL(eventQC, processQC<FilteredMyCollisions_Qvec>, processEventQC_Cent_Qvec, "event QC + q vector", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<eventQC>(cfgc, TaskName{"event-qc"})};
}
