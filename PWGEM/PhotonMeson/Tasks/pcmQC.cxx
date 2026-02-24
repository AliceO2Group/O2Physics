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

/// \file MaterialBudget.cxx
/// \brief This code runs loop over v0 photons for PCM QC.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch

#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <MathUtils/Utils.h>

#include <TH1.h>

#include <cmath>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents_004, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyV0PhotonsML = soa::Join<aod::V0PhotonsKF, aod::V0PhotonsPhiVPsi, aod::V0KFEMEventIds>;
using MyV0PhotonML = MyV0PhotonsML::iterator;

struct PCMQC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_pt_v0{"cfg_max_pt_v0", 1e+10, "max pT for v0 photons at PV"};
    Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", +0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 1.5, "max distance btween 2 legs"};
    Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};
    Configurable<bool> cfg_dEdx_postcalibration{"cfg_dEdx_postcalibration", false, "flag to enable dEdx post calibration"};
    // for ML cuts
    Configurable<bool> cfg_apply_ml_cuts{"cfg_apply_ml", false, "flag to apply ML cut"};
    Configurable<bool> cfg_use_2d_binning{"cfg_use_2d_binning", false, "flag to use 2D binning (pT, cent)"};
    Configurable<bool> cfg_load_ml_models_from_ccdb{"cfg_load_ml_models_from_ccdb", true, "flag to load ML models from CCDB"};
    Configurable<int> cfg_timestamp_ccdb{"cfg_timestamp_ccdb", -1, "timestamp for CCDB"};
    Configurable<int> cfg_nclasses_ml{"cfg_nclasses_ml", static_cast<int>(o2::analysis::em_cuts_ml::NCutScores), "number of classes for ML"};
    Configurable<std::vector<int>> cfg_cut_dir_ml{"cfg_cut_dir_ml", std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}, "cut direction for ML"};
    Configurable<std::vector<std::string>> cfg_input_feature_names{"cfg_input_feature_names", std::vector<std::string>{"feature1", "feature2"}, "input feature names for ML models"};
    Configurable<std::vector<std::string>> cfg_model_paths_ccdb{"cfg_model_paths_ccdb", std::vector<std::string>{"path_ccdb/BDT_PCM/"}, "CCDB paths for ML models"};
    Configurable<std::vector<std::string>> cfg_onnx_file_names{"cfg_onnx_file_names", std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}, "ONNX file names for ML models"};
    Configurable<std::vector<std::string>> cfg_labels_bins_ml{"cfg_labels_bins_ml", std::vector<std::string>{"bin 0", "bin 1"}, "Labels for bins"};
    Configurable<std::vector<std::string>> cfg_labels_cut_scores_ml{"cfg_labels_cut_scores_ml", std::vector<std::string>{o2::analysis::em_cuts_ml::labelsCutScore}, "Labels for cut scores"};
    Configurable<std::vector<double>> cfg_bins_pt_ml{"cfg_bins_pt_ml", std::vector<double>{0.0, +1e+10}, "pT bin limits for ML application"};
    Configurable<std::vector<double>> cfg_bins_cent_ml{"cfg_bins_cent_ml", std::vector<double>{o2::analysis::em_cuts_ml::vecBinsCent}, "centrality bins for ML"};
    Configurable<std::vector<double>> cfg_cuts_ml_flat{"cfg_cuts_ml_flat", {0.5}, "Flattened ML cuts: [bin0_score0, bin0_score1, ..., binN_scoreM]"};
  } pcmcuts;

  o2::ccdb::CcdbApi ccdbApi;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    addhistograms();
    DefineEMEventCut();
    DefinePCMCut();

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    fV0PhotonCut.SetD_Bz(d_bz);
    mRunNumber = collision.runNumber();
  }

  void addhistograms()
  {
    // event info
    auto hCollisionCounter = fRegistry.add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{10, 0.5, 10.5}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(10, "accepted");

    fRegistry.add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{300, 0, 6000}, {300, 0, 6000}}, false);
    fRegistry.add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0MvsMultNTracksPV", "hCentFT0MvsMultNTracksPV;centrality FT0M (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
    fRegistry.add("Event/before/hMultFT0MvsMultNTracksPV", "hMultFT0MvsMultNTracksPV;mult. FT0M;N_{track} to PV", kTH2F, {{600, 0, 6000}, {600, 0, 6000}}, false);
    fRegistry.addClone("Event/before/", "Event/after/");

    // v0 info
    fRegistry.add("V0/hPt", "pT;p_{T,#gamma} (GeV/c)", kTH1F, {{2000, 0.0f, 20}}, false);
    fRegistry.add("V0/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, o2::constants::math::TwoPI}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("V0/hXY", "conversion point in XY;V_{x} (cm);V_{y} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, -100.0f, 100.0f}}, false);
    fRegistry.add("V0/hRZ", "conversion point in RZ;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100, 100}, {200, 0.0f, 100.0f}}, false);
    fRegistry.add("V0/hCosPA", "V0CosPA;cosine pointing angle in 3D", kTH1F, {{100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/hCosPAXY", "V0CosPA;cosine pointing angle in XY", kTH1F, {{100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/hCosPARZ", "V0CosPA;cosine pointing angle in RZ", kTH1F, {{100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/hPCA", "distance between 2 legs;PCA (cm)", kTH1F, {{500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}, false);
    fRegistry.add("V0/hDCAz_Pt", "DCA_{z} to PV vs. p_{T};DCA_{z} (cm);p_{T} (GeV/c)", kTH2F, {{200, -5.f, +5.f}, {2000, 0.0f, 20}}, false);
    fRegistry.add("V0/hAPplot", "AP plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}}, false);
    fRegistry.add("V0/hMassGamma", "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.0f, 0.1f}}, false);
    fRegistry.add("V0/hKFChi2vsM", "KF chi2 vs. m_{ee};m_{ee} (GeV/c^{2});KF chi2/NDF", kTH2F, {{100, 0.0f, 0.1f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsR", "KF chi2 vs. conversion point in XY;R_{xy} (cm);KF chi2/NDF", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsX", "KF chi2 vs. conversion point in X;X (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsY", "KF chi2 vs. conversion point in Y;Y (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsZ", "KF chi2 vs. conversion point in Z;Z (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hsConvPoint", "photon conversion point;r_{xy} (cm);#varphi (rad.);#eta;", kTHnSparseF, {{100, 0.0f, 100}, {90, 0, o2::constants::math::TwoPI}, {80, -2, +2}}, false);
    fRegistry.add("V0/hNgamma", "Number of #gamma candidates per collision", kTH1F, {{101, -0.5f, 100.5f}});

    if (pcmcuts.cfg_apply_ml_cuts) {
      if (pcmcuts.cfg_nclasses_ml == 2) {
        fRegistry.add("V0/hBDTBackgroundScoreVsPt", "BDT background score vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hBDTSignalScoreVsPt", "BDT signal score vs pT; pT (GeV/c); BDT signal score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hPhiVPsi", "#varphi vs. #psi angle;#psi (rad.); #varphi (rad.)", kTH2F, {{200, -o2::constants::math::PI, o2::constants::math::PI}, {200, 0, o2::constants::math::TwoPI}}, false);
      } else if (pcmcuts.cfg_nclasses_ml == 3) {
        fRegistry.add("V0/hBDTBackgroundScoreVsPt", "BDT background score vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hBDTPrimaryPhotonScoreVsPt", "BDT primary photon score vs pT; pT (GeV/c); BDT primary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hBDTSecondaryPhotonScoreVsPt", "BDT secondary photon score vs pT; pT (GeV/c); BDT secondary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hPhiVPsi", "#varphi vs. #psi angle;#psi (rad.); #varphi (rad.)", kTH2F, {{200, -o2::constants::math::PI, o2::constants::math::PI}, {200, 0, o2::constants::math::TwoPI}}, false);
      } else {
        fRegistry.add("V0/hBDTScoreVsPt", "BDT score vs pT; pT (GeV/c); BDT score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        fRegistry.add("V0/hPhiVPsi", "#varphi vs. #psi angle;#psi (rad.); #varphi (rad.)", kTH2F, {{200, -o2::constants::math::PI, o2::constants::math::PI}, {200, 0, o2::constants::math::TwoPI}}, false);
      }
    }

    // v0leg info
    fRegistry.add("V0Leg/hPt", "pT;p_{T,e} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("V0Leg/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -50, 50}}, false);
    fRegistry.add("V0Leg/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, o2::constants::math::TwoPI}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("V0Leg/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -50.0f, 50.0f}, {200, -50.0f, 50.0f}}, false);
    fRegistry.add("V0Leg/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("V0Leg/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    fRegistry.add("V0Leg/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("V0Leg/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("V0Leg/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {160, 0, 16}}, false);
    if (pcmcuts.cfg_dEdx_postcalibration) {
      fRegistry.add("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Pos", "momentum of pos leg vs. conversion point of V0 vs. TPC n sigma pos vs. eta of pos leg; p (GeV/c); r_{xy} (cm); n #sigma_{e}^{TPC}; #eta", kTHnSparseF, {{200, 0, 20}, {100, 0, 100}, {500, -5, 5}, {200, -1, +1}}, false);
      fRegistry.add("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Ele", "momentum of neg leg vs. conversion point of V0 vs. TPC n sigma el vs. eta of neg leg; p (GeV/c); r_{xy} (cm); n #sigma_{e}^{TPC}; #eta", kTHnSparseF, {{200, 0, 20}, {100, 0, 100}, {500, -5, 5}, {200, -1, +1}}, false);
    }
    // fRegistry.add("V0Leg/hXY", "X vs. Y;X (cm);Y (cm)", kTH2F, {{100, 0, 100}, {80, -20, 20}}, false);
    // fRegistry.add("V0Leg/hZX", "Z vs. X;Z (cm);X (cm)", kTH2F, {{200, -100, 100}, {100, 0, 100}}, false);
    // fRegistry.add("V0Leg/hZY", "Z vs. Y;Z (cm);Y (cm)", kTH2F, {{200, -100, 100}, {80, -20, 20}}, false);
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, pcmcuts.cfg_max_pt_v0);
    fV0PhotonCut.SetV0EtaRange(pcmcuts.cfg_min_eta_v0, pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfg_max_chi2kf);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfg_max_frac_shared_clusters_tpc);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfg_disable_itsonly_track);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfg_disable_tpconly_track);
    fV0PhotonCut.SetRequireITSTPC(pcmcuts.cfg_require_v0_with_itstpc);
    fV0PhotonCut.SetRequireITSonly(pcmcuts.cfg_require_v0_with_itsonly);
    fV0PhotonCut.SetRequireTPConly(pcmcuts.cfg_require_v0_with_tpconly);

    // for ML
    fV0PhotonCut.SetApplyMlCuts(pcmcuts.cfg_apply_ml_cuts);
    fV0PhotonCut.SetUse2DBinning(pcmcuts.cfg_use_2d_binning);
    fV0PhotonCut.SetLoadMlModelsFromCCDB(pcmcuts.cfg_load_ml_models_from_ccdb);
    fV0PhotonCut.SetNClassesMl(pcmcuts.cfg_nclasses_ml);
    fV0PhotonCut.SetMlTimestampCCDB(pcmcuts.cfg_timestamp_ccdb);
    fV0PhotonCut.SetCcdbUrl(ccdburl);
    CentType mCentralityTypeMlEnum;
    mCentralityTypeMlEnum = static_cast<CentType>(cfgCentEstimator.value);
    fV0PhotonCut.SetCentralityTypeMl(mCentralityTypeMlEnum);
    fV0PhotonCut.SetCutDirMl(pcmcuts.cfg_cut_dir_ml);
    fV0PhotonCut.SetMlModelPathsCCDB(pcmcuts.cfg_model_paths_ccdb);
    fV0PhotonCut.SetMlOnnxFileNames(pcmcuts.cfg_onnx_file_names);
    fV0PhotonCut.SetBinsPtMl(pcmcuts.cfg_bins_pt_ml);
    fV0PhotonCut.SetBinsCentMl(pcmcuts.cfg_bins_cent_ml);
    fV0PhotonCut.SetCutsMl(pcmcuts.cfg_cuts_ml_flat);
    fV0PhotonCut.SetNamesInputFeatures(pcmcuts.cfg_input_feature_names);
    fV0PhotonCut.SetLabelsBinsMl(pcmcuts.cfg_labels_bins_ml);
    fV0PhotonCut.SetLabelsCutScoresMl(pcmcuts.cfg_labels_cut_scores_ml);
    fV0PhotonCut.SetD_Bz(0.0f); // dummy value -> only for psi_pair calculation

    if (pcmcuts.cfg_apply_ml_cuts) {
      fV0PhotonCut.initV0MlModels(ccdbApi);
    }
  }

  template <const int ev_id, typename TCollision>
  void fillEventInfo(TCollision const& collision, const float /*weight*/ = 1.f)
  {
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
    }
    if (collision.sel8()) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
    }
    if (std::fabs(collision.posZ()) < 10.0) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
    }
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0MvsMultNTracksPV"), collision.centFT0M(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0MvsMultNTracksPV"), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV());
  }

  template <typename TV0>
  void fillV0Info(TV0 const& v0)
  {
    fRegistry.fill(HIST("V0/hPt"), v0.pt());
    fRegistry.fill(HIST("V0/hEtaPhi"), v0.phi(), v0.eta());
    fRegistry.fill(HIST("V0/hXY"), v0.vx(), v0.vy());
    fRegistry.fill(HIST("V0/hRZ"), v0.vz(), v0.v0radius());
    fRegistry.fill(HIST("V0/hCosPA"), v0.cospa());
    fRegistry.fill(HIST("V0/hCosPAXY"), v0.cospaXY());
    fRegistry.fill(HIST("V0/hCosPARZ"), v0.cospaRZ());
    fRegistry.fill(HIST("V0/hPCA"), v0.pca());
    fRegistry.fill(HIST("V0/hDCAxyz"), v0.dcaXYtopv(), v0.dcaZtopv());
    fRegistry.fill(HIST("V0/hDCAz_Pt"), v0.dcaZtopv(), v0.pt());
    fRegistry.fill(HIST("V0/hAPplot"), v0.alpha(), v0.qtarm());
    fRegistry.fill(HIST("V0/hMassGamma"), v0.v0radius(), v0.mGamma());
    fRegistry.fill(HIST("V0/hKFChi2vsM"), v0.mGamma(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsR"), v0.v0radius(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsX"), v0.vx(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsY"), v0.vy(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsZ"), v0.vz(), v0.chiSquareNDF());

    float phi_cp = std::atan2(v0.vy(), v0.vx());
    o2::math_utils::bringTo02Pi(phi_cp);
    float eta_cp = std::atanh(v0.vz() / std::sqrt(std::pow(v0.vx(), 2) + std::pow(v0.vy(), 2) + std::pow(v0.vz(), 2)));
    fRegistry.fill(HIST("V0/hsConvPoint"), v0.v0radius(), phi_cp, eta_cp);

    // BDT response histogram can be filled here when apply BDT is true
    if (pcmcuts.cfg_apply_ml_cuts) {
      const std::span<const float>& bdtValue = fV0PhotonCut.getBDTValue();
      float psipair = 999.f;
      float phiv = 999.f;
      if constexpr (requires { v0.psipair(); v0.phiv(); }) {
        psipair = v0.psipair();
        phiv = v0.phiv();
      }
      fRegistry.fill(HIST("V0/hPhiVPsi"), psipair, phiv);
      if (pcmcuts.cfg_nclasses_ml == 2 && bdtValue.size() == 2) {
        fRegistry.fill(HIST("V0/hBDTBackgroundScoreVsPt"), v0.pt(), bdtValue[0]);
        fRegistry.fill(HIST("V0/hBDTSignalScoreVsPt"), v0.pt(), bdtValue[1]);
      } else if (pcmcuts.cfg_nclasses_ml == 3 && bdtValue.size() == 3) {
        fRegistry.fill(HIST("V0/hBDTBackgroundScoreVsPt"), v0.pt(), bdtValue[0]);
        fRegistry.fill(HIST("V0/hBDTPrimaryPhotonScoreVsPt"), v0.pt(), bdtValue[1]);
        fRegistry.fill(HIST("V0/hBDTSecondaryPhotonScoreVsPt"), v0.pt(), bdtValue[2]);
      } else if (bdtValue.size() == 1) {
        fRegistry.fill(HIST("V0/hBDTCutVsPt"), v0.pt(), bdtValue[0]);
      }
    }
  }

  template <typename TLeg>
  void fillV0LegInfo(TLeg const& leg)
  {
    fRegistry.fill(HIST("V0Leg/hPt"), leg.pt());
    fRegistry.fill(HIST("V0Leg/hQoverPt"), leg.sign() / leg.pt());
    fRegistry.fill(HIST("V0Leg/hEtaPhi"), leg.phi(), leg.eta());
    fRegistry.fill(HIST("V0Leg/hDCAxyz"), leg.dcaXY(), leg.dcaZ());
    fRegistry.fill(HIST("V0Leg/hNclsITS"), leg.itsNCls());
    fRegistry.fill(HIST("V0Leg/hNclsTPC"), leg.tpcNClsFound());
    fRegistry.fill(HIST("V0Leg/hNcrTPC"), leg.tpcNClsCrossedRows());
    fRegistry.fill(HIST("V0Leg/hTPCNcr2Nf"), leg.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("V0Leg/hTPCNcls2Nf"), leg.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("V0Leg/hTPCNclsShared"), leg.pt(), leg.tpcFractionSharedCls());
    fRegistry.fill(HIST("V0Leg/hChi2TPC"), leg.tpcChi2NCl());
    fRegistry.fill(HIST("V0Leg/hChi2ITS"), leg.itsChi2NCl());
    fRegistry.fill(HIST("V0Leg/hITSClusterMap"), leg.itsClusterMap());
    if (leg.hasITS()) {
      fRegistry.fill(HIST("V0Leg/hMeanClusterSizeITS"), leg.p(), leg.meanClusterSizeITS() * std::cos(std::atan(leg.tgl())));
    }
    fRegistry.fill(HIST("V0Leg/hTPCdEdx"), leg.tpcInnerParam(), leg.tpcSignal());
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaPi"), leg.tpcInnerParam(), leg.tpcNSigmaPi());
    // fRegistry.fill(HIST("V0Leg/hXY"), leg.x(), leg.y());
    // fRegistry.fill(HIST("V0Leg/hZX"), leg.z(), leg.x());
    // fRegistry.fill(HIST("V0Leg/hZY"), leg.z(), leg.y());
  }

  o2::framework::SliceCache v0cache;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  template <typename TV0Photon>
  void process(FilteredMyCollisions const& collisions, TV0Photon const& v0photons, aod::V0Legs const&)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      fillEventInfo<0>(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      fillEventInfo<1>(collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      fV0PhotonCut.SetCentrality(centralities[cfgCentEstimator]);
      int nv0 = 0;
      auto v0photons_coll = v0photons.sliceByCached(aod::v0photonkf::emeventId, collision.globalIndex(), v0cache);
      for (const auto& v0 : v0photons_coll) {
        auto pos = v0.template posTrack_as<aod::V0Legs>();
        auto ele = v0.template negTrack_as<aod::V0Legs>();

        if (!fV0PhotonCut.IsSelected<decltype(v0), aod::V0Legs>(v0)) {
          continue;
        }
        fillV0Info(v0);
        for (auto& leg : {pos, ele}) {
          fillV0LegInfo(leg);
        }
        if (pcmcuts.cfg_dEdx_postcalibration) {
          fRegistry.fill(HIST("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Pos"), pos.p(), v0.v0radius(), pos.tpcNSigmaEl(), pos.eta());
          fRegistry.fill(HIST("V0Leg/hPvsConvPointvsTPCNsigmaElvsEta_Ele"), ele.p(), v0.v0radius(), ele.tpcNSigmaEl(), ele.eta());
        }
        nv0++;
      } // end of v0 loop
      fRegistry.fill(HIST("V0/hNgamma"), nv0);
    } // end of collision loop
  }
  void processQC(FilteredMyCollisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const& v0legs)
  {
    process(collisions, v0photons, v0legs);
  } // end of process

  void processQCML(FilteredMyCollisions const& collisions, MyV0PhotonsML const& v0photonsML, aod::V0Legs const& v0legs)
  {
    process(collisions, v0photonsML, v0legs);
  } // end of ML process

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(PCMQC, processQC, "run PCM QC", true);
  PROCESS_SWITCH(PCMQC, processQCML, "run PCM QC with ML", false);
  PROCESS_SWITCH(PCMQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQC>(cfgc, TaskName{"pcm-qc"})};
}
