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
//
// Analysis task to produce resolution mapfor electrons/muons in dilepton analysis
//    Please write to: daiki.sekihata@cern.ch

#include <map>
#include <string>
#include <utility>
#include <set>
#include <vector>
#include <array>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/RCTSelectionFlags.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"

#include "DetectorsBase/Propagator.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "Common/Core/fwdtrackUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::fwdtrackutils;

struct CreateResolutionMap {
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<bool> cfg_require_true_mc_collision_association{"cfg_require_true_mc_collision_association", false, "flag to require true mc collision association"};
  Configurable<bool> cfg_reject_fake_match_its_tpc{"cfg_reject_fake_match_its_tpc", false, "flag to reject fake match between ITS-TPC"};
  // Configurable<bool> cfg_reject_fake_match_its_tpc_tof{"cfg_reject_fake_match_its_tpc_tof", false, "flag to reject fake match between ITS-TPC-TOF"};
  Configurable<bool> cfg_reject_fake_match_mft_mch{"cfg_reject_fake_match_mft_mch", false, "flag to reject fake match between MFT-MCH"};

  ConfigurableAxis ConfPtGenBins{"ConfPtGenBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00}, "gen. pT bins for output histograms"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0, 10, 30, 50, 110}, "centrality (%) bins for output histograms"};

  ConfigurableAxis ConfEtaCBGenBins{"ConfEtaCBGenBins", {30, -1.5, +1.5}, "gen. eta bins at midrapidity for output histograms"};
  ConfigurableAxis ConfEtaFWDGenBins{"ConfEtaFWDGenBins", {40, -5.5, -1.5}, "gen. eta bins at forward rapidity for output histograms"};
  ConfigurableAxis ConfPhiGenBins{"ConfPhiGenBins", {72, 0, 2.f * M_PI}, "gen. eta bins at forward rapidity for output histograms"};

  ConfigurableAxis ConfRelDeltaPtBins{"ConfRelDeltaPtBins", {200, -1.f, +1.f}, "rel. dpt for output histograms"};
  ConfigurableAxis ConfDeltaEtaBins{"ConfDeltaEtaBins", {200, -0.2f, +0.2f}, "deta bins for output histograms"};
  ConfigurableAxis ConfDeltaPhiBins{"ConfDeltaPhiBins", {200, -0.2f, +0.2f}, "dphi bins for output histograms"};

  Configurable<bool> cfgFillTHnSparse{"cfgFillTHnSparse", true, "fill THnSparse for output"};
  Configurable<bool> cfgFillTH2{"cfgFillTH2", false, "fill TH2 for output"};

  Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
  Configurable<std::string> cfgRCTLabelCB{"cfgRCTLabelCB", "CBT_hadronPID", "select 1 [CBT, CBT_hadron] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<std::string> cfgRCTLabelFWDSA{"cfgRCTLabelFWDSA", "CBT_muon", "select 1 [CBT_muon] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<std::string> cfgRCTLabelFWDGL{"cfgRCTLabelFWDGL", "CBT_muon_glo", "select 1 [CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
  Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    // Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    // Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    // Configurable<bool> cfgRequirekNoCollInRofStandard{"cfgRequirekNoCollInRofStandard", false, "require no other collisions in this Readout Frame with per-collision multiplicity above threshold"};
    // Configurable<bool> cfgRequirekNoCollInRofStrict{"cfgRequirekNoCollInRofStrict", false, "require no other collisions in this Readout Frame"};
    // Configurable<bool> cfgRequirekNoHighMultCollInPrevRof{"cfgRequirekNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
    // Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "number of inactive chips on ITS layer 3 are below threshold "};
    // Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "number of inactive chips on ITS layers 0-3 are below threshold "};
    // Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "number of inactive chips on all ITS layers are below threshold "};
  } eventcuts;

  struct : ConfigurableGroup {
    std::string prefix = "electroncut_group";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -1.5, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +1.5, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 80, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_tpc_cr_findable_ratio{"cfg_min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
  } electroncuts;

  struct : ConfigurableGroup {
    std::string prefix = "muoncut_group";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
    Configurable<float> cfg_min_eta_track_sa{"cfg_min_eta_track_sa", -5.5, "min eta for standalone muon track"};
    Configurable<float> cfg_max_eta_track_sa{"cfg_max_eta_track_sa", -1.5, "max eta for standalone muon track"};
    Configurable<float> cfg_min_eta_track_gl{"cfg_min_eta_track_gl", -5.5, "min eta for global muon track"};
    Configurable<float> cfg_max_eta_track_gl{"cfg_max_eta_track_gl", -1.5, "max eta for global muon track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2_sa{"cfg_max_chi2_sa", 1e+10, "max chi2 for standalone muon track"};
    Configurable<float> cfg_max_chi2_gl{"cfg_max_chi2_gl", 40, "max chi2 for standalone muon track"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 40, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy_gl{"cfg_max_dcaxy_gl", 0.1, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs_sa{"cfg_min_rabs_sa", 17.6, "min Radius at the absorber end for standalone muon track"};
    Configurable<float> cfg_max_rabs_sa{"cfg_max_rabs_sa", 89.5, "max Radius at the absorber end for standalone muon track"};
    Configurable<float> cfg_min_rabs_gl{"cfg_min_rabs_gl", 27.6, "min Radius at the absorber end for global muon track"};
    Configurable<float> cfg_max_rabs_gl{"cfg_max_rabs_gl", 89.5, "max Radius at the absorber end for global muon track"};
    Configurable<float> cfg_mid_rabs{"cfg_mid_rabs", 26.5, "middle R at absorber end for pDCA cut"};
    Configurable<float> cfg_max_pdca_forLargeR{"cfg_max_pdca_forLargeR", 324.f, "max. pDCA for large R at absorber end"};
    Configurable<float> cfg_max_pdca_forSmallR{"cfg_max_pdca_forSmallR", 594.f, "max. pDCA for small R at absorber end"};
    Configurable<int> cfg_min_nclusters_MFT{"cfg_min_nclusters_MFT", 5, "min nclusters MFT"};
    Configurable<int> cfg_min_nclusters_MCH{"min_min_nclusters_MCH", 5, "min nclusters MCH"};
    Configurable<float> cfg_max_reldpt{"cfg_max_reldpt", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_deta{"cfg_max_deta", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_dphi{"cfg_max_dphi", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
    Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
  } muoncuts;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::aod::rctsel::RCTFlagsChecker rctCheckerCB;
  o2::aod::rctsel::RCTFlagsChecker rctCheckerFWDSA;
  o2::aod::rctsel::RCTFlagsChecker rctCheckerFWDGL;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float d_bz;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  void init(o2::framework::InitContext&)
  {
    if (doprocessElectronSA && doprocessElectronTTCA) {
      LOGF(fatal, "Cannot enable processElectronSA and processElectronTTCA at the same time. Please choose one.");
    }

    if (doprocessMuonSA && doprocessMuonTTCA) {
      LOGF(fatal, "Cannot enable processMuonSA and processMuonTTCA at the same time. Please choose one.");
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    rctCheckerCB.init(cfgRCTLabelCB.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);
    rctCheckerFWDSA.init(cfgRCTLabelFWDSA.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);
    rctCheckerFWDGL.init(cfgRCTLabelFWDGL.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);

    mRunNumber = 0;
    d_bz = 0;

    const AxisSpec axis_cent{ConfCentBins, "centrality (%)"};
    const AxisSpec axis_pt_gen{ConfPtGenBins, "p_{T,l}^{gen} (GeV/c)"};
    const AxisSpec axis_eta_cb_gen{ConfEtaCBGenBins, "#eta_{l}^{gen}"};
    const AxisSpec axis_eta_fwd_gen{ConfEtaFWDGenBins, "#eta_{l}^{gen}"};
    const AxisSpec axis_phi_gen{ConfPhiGenBins, "#varphi_{l}^{gen} (rad.)"};
    const AxisSpec axis_dpt{ConfRelDeltaPtBins, "(p_{T,l}^{gen} - p_{T,l}^{rec})/p_{T,l}^{gen}"};
    const AxisSpec axis_deta{ConfDeltaEtaBins, "#eta_{l}^{gen} - #eta_{l}^{rec}"};
    const AxisSpec axis_dphi{ConfDeltaPhiBins, "#varphi_{l}^{gen} - #varphi_{l}^{rec} (rad.)"};
    const AxisSpec axis_charge_gen{3, -1.5, +1.5, "true sign"};

    // registry.add("Event/Electron/hImpPar_Centrality", "true imapact parameter vs. estimated centrality;impact parameter (fm);centrality (%)", kTH2F, {{200, 0, 20}, {110, 0, 110}}, true);
    // registry.add("Event/Electron/hImpPar_Centrality", "true imapact parameter vs. estimated centrality;impact parameter (fm);centrality (%)", kTH2F, {{200, 0, 20}, {110, 0, 110}}, true);
    if (doprocessGen) {
      registry.add("Event/hGenID", "generator ID;generator ID;Number of mc collisions", kTH1F, {{7, -1.5, 5.5}}, true);
    }
    if (cfgFillTH2) {
      registry.add("Electron/hPt", "rec. p_{T,l};p_{T,l} (GeV/c)", kTH1F, {{1000, 0, 10}}, false);
      registry.add("Electron/hEtaPhi", "rec. #eta vs. #varphi;#varphi_{l} (rad.);#eta_{l}", kTH2F, {{90, 0, 2 * M_PI}, {100, -5, +5}}, false);
      registry.add("Electron/Ptgen_RelDeltaPt", "resolution", kTH2F, {{axis_pt_gen}, {axis_dpt}}, true);
      registry.add("Electron/Ptgen_DeltaEta", "resolution", kTH2F, {{axis_pt_gen}, {axis_deta}}, true);
      registry.add("Electron/Ptgen_DeltaPhi_Pos", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
      registry.add("Electron/Ptgen_DeltaPhi_Neg", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
      registry.addClone("Electron/", "StandaloneMuon/");
      registry.addClone("Electron/", "GlobalMuon/");
    }

    if (cfgFillTHnSparse) {
      registry.add("Electron/hs_reso", "8D resolution", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_cb_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
      registry.add("StandaloneMuon/hs_reso", "8D resolution", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
      registry.add("GlobalMuon/hs_reso", "8D resolution", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      mRunNumber = bc.runNumber();

      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
      o2::mch::TrackExtrap::setField();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    }
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());

      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();

    // std::map<string, string> metadata;
    // auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    // auto ts = soreor.first;
    // auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    // o2::base::Propagator::initFieldFromGRP(grpmag);

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
  }

  template <typename TCollision>
  bool isSelectedEvent(TCollision const& collision)
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

    if (eventcuts.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (eventcuts.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (!(eventcuts.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgTrackOccupancyMax)) {
      return false;
    }

    if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      return false;
    }

    // if (eventcuts.cfgRequireNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
    //   return false;
    // }

    // if (eventcuts.cfgRequireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
    //   return false;
    // }

    // if (eventcuts.cfgRequirekNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
    //   return false;
    // }

    // if (eventcuts.cfgRequirekNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
    //   return false;
    // }

    // if (eventcuts.cfgRequirekNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
    //   return false;
    // }

    // if (eventcuts.cfgRequireGoodITSLayer3 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
    //   return false;
    // }

    // if (eventcuts.cfgRequireGoodITSLayer0123 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
    //   return false;
    // }

    // if (eventcuts.cfgRequireGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
    //   return false;
    // }

    return true;
  }

  std::pair<int8_t, std::set<uint8_t>> itsRequirement_ibany = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.
  std::pair<int8_t, std::set<uint8_t>> itsRequirement_ib1st = {1, {0}};       // first hit on ITS ib layers.

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.tpcChi2NCl() > electroncuts.cfg_max_chi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > electroncuts.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < electroncuts.cfg_min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < electroncuts.cfg_min_ncluster_itsib) {
      return false;
    }

    auto hits = std::count_if(itsRequirement_ibany.second.begin(), itsRequirement_ibany.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement_ibany.first) {
      return false;
    }
    if (electroncuts.cfg_require_itsib_1st) {
      auto hit_ib1st = std::count_if(itsRequirement_ib1st.second.begin(), itsRequirement_ib1st.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      if (hit_ib1st < itsRequirement_ib1st.first) {
        return false;
      }
    }

    if (track.tpcNClsFound() < electroncuts.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < electroncuts.cfg_min_ncrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < electroncuts.cfg_min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcFractionSharedCls() > electroncuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isSelectedTrackKine(TTrack const& track, const float pt, const float eta, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.tpcChi2NCl() > electroncuts.cfg_max_chi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > electroncuts.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < electroncuts.cfg_min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < electroncuts.cfg_min_ncluster_itsib) {
      return false;
    }

    auto hits = std::count_if(itsRequirement_ibany.second.begin(), itsRequirement_ibany.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement_ibany.first) {
      return false;
    }
    if (electroncuts.cfg_require_itsib_1st) {
      auto hit_ib1st = std::count_if(itsRequirement_ib1st.second.begin(), itsRequirement_ib1st.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      if (hit_ib1st < itsRequirement_ib1st.first) {
        return false;
      }
    }

    if (track.tpcNClsFound() < electroncuts.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < electroncuts.cfg_min_ncrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < electroncuts.cfg_min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcFractionSharedCls() > electroncuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    if (std::fabs(dcaXY) > electroncuts.cfg_max_dcaxy || std::fabs(dcaZ) > electroncuts.cfg_max_dcaz) {
      return false;
    }

    if (pt < electroncuts.cfg_min_pt_track || std::fabs(eta) > electroncuts.cfg_max_eta_track) {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TMuon>
  void fillMuon(TCollision const& collision, TMuon const& muon, const float centrality)
  {
    auto mcparticle = muon.template mcParticle_as<aod::McParticles>();
    if (std::abs(mcparticle.pdgCode()) != 13 || !(mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator())) {
      return;
    }
    if (cfg_require_true_mc_collision_association && mcparticle.mcCollisionId() != collision.mcCollisionId()) {
      return;
    }
    if (muon.trackType() != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) && muon.trackType() != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      return;
    }

    if (muon.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return;
    }
    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(muon, collision, propagationPoint::kToVertex);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = propagateMuon(muon, collision, propagationPoint::kToDCA);

    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    float dcaX = propmuonAtDCA.getX() - collision.posX();
    float dcaY = propmuonAtDCA.getY() - collision.posY();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

    float rAtAbsorberEnd = muon.rAtAbsorberEnd(); // this works only for GlobalMuonTrack
    float pDCA = muon.p() * dcaXY;
    int nClustersMFT = 0;
    float ptMatchedMCHMID = propmuonAtPV.getPt();
    float etaMatchedMCHMID = propmuonAtPV.getEta();
    float phiMatchedMCHMID = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);

    if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      // mcparticle for global MFT-MCH-MID is identical to mcparticle of MCH-MID track. If not, mismatch.
      const auto& mchtrack = muon.template matchMCHTrack_as<MyFwdTracks>(); // MCH-MID
      const auto& mfttrack = muon.template matchMFTTrack_as<MyMFTTracks>(); // MFTsa
      if (!mchtrack.has_mcParticle() || !mfttrack.has_mcParticle()) {
        return;
      }
      auto mcparticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>();
      auto mcparticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
      if (mcparticle.globalIndex() != mcparticle_MCHMID.globalIndex()) { // this should not happen. this is only for protection.
        return;
      }
      if (cfg_reject_fake_match_mft_mch && mcparticle.globalIndex() != mcparticle_MFT.globalIndex()) { // evaluate mismatch
        return;
      }

      o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, collision, propagationPoint::kToVertex);
      ptMatchedMCHMID = propmuonAtPV_Matched.getPt();
      etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
      phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
      o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
      o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, collision, propagationPoint::kToDCA);
      float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
      float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
      float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
      pDCA = mchtrack.p() * dcaXY_Matched;
      nClustersMFT = mfttrack.nClusters();

      if (nClustersMFT < muoncuts.cfg_min_ncluster_mft) {
        return;
      }
      if (muon.chi2MatchMCHMFT() > muoncuts.cfg_max_matching_chi2_mftmch) {
        return;
      }
      if (muoncuts.refitGlobalMuon) {
        eta = mfttrack.eta();
        phi = mfttrack.phi();
        o2::math_utils::bringTo02Pi(phi);
        pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
      }

      float dpt = (ptMatchedMCHMID - pt) / pt;
      float deta = etaMatchedMCHMID - eta;
      float dphi = phiMatchedMCHMID - phi;
      o2::math_utils::bringToPMPi(dphi);
      if (std::sqrt(std::pow(deta / muoncuts.cfg_max_deta, 2) + std::pow(dphi / muoncuts.cfg_max_dphi, 2)) > 1.f || std::fabs(dpt) > muoncuts.cfg_max_reldpt) {
        return;
      }

    } else if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = propagateMuon(muon, collision, propagationPoint::kToRabs); // this is necessary only for MuonStandaloneTrack
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    } else {
      return;
    }

    if (muon.nClusters() < muoncuts.cfg_min_ncluster_mch) {
      return;
    }

    if (!isSelectedMuon(pt, eta, rAtAbsorberEnd, pDCA, muon.chi2(), muon.trackType(), dcaXY)) {
      return;
    }

    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      if (cfgRequireGoodRCT && !rctCheckerFWDSA.checkTable(collision)) {
        return;
      }
      if (cfgFillTHnSparse) {
        registry.fill(HIST("StandaloneMuon/hs_reso"), centrality, mcparticle.pt(), mcparticle.eta(), mcparticle.phi(), -mcparticle.pdgCode() / 13, (mcparticle.pt() - pt) / mcparticle.pt(), mcparticle.eta() - eta, mcparticle.phi() - phi);
      }

      if (cfgFillTH2) {
        registry.fill(HIST("StandaloneMuon/hPt"), pt);
        registry.fill(HIST("StandaloneMuon/hEtaPhi"), phi, eta);
        registry.fill(HIST("StandaloneMuon/Ptgen_RelDeltaPt"), mcparticle.pt(), (mcparticle.pt() - pt) / mcparticle.pt());
        registry.fill(HIST("StandaloneMuon/Ptgen_DeltaEta"), mcparticle.pt(), mcparticle.eta() - eta);
        if (mcparticle.pdgCode() == -13) { // positive muon
          registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Pos"), mcparticle.pt(), mcparticle.phi() - phi);
        } else if (mcparticle.pdgCode() == 13) { // negative muon
          registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Neg"), mcparticle.pt(), mcparticle.phi() - phi);
        }
      }
    } else if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      if (cfgRequireGoodRCT && !rctCheckerFWDGL.checkTable(collision)) {
        return;
      }
      if (cfgFillTHnSparse) {
        registry.fill(HIST("GlobalMuon/hs_reso"), centrality, mcparticle.pt(), mcparticle.eta(), mcparticle.phi(), -mcparticle.pdgCode() / 13, (mcparticle.pt() - pt) / mcparticle.pt(), mcparticle.eta() - eta, mcparticle.phi() - phi);
      }
      if (cfgFillTH2) {
        registry.fill(HIST("GlobalMuon/hPt"), pt);
        registry.fill(HIST("GlobalMuon/hEtaPhi"), phi, eta);
        registry.fill(HIST("GlobalMuon/Ptgen_RelDeltaPt"), mcparticle.pt(), (mcparticle.pt() - pt) / mcparticle.pt());
        registry.fill(HIST("GlobalMuon/Ptgen_DeltaEta"), mcparticle.pt(), mcparticle.eta() - eta);
        if (mcparticle.pdgCode() == -13) { // positive muon
          registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Pos"), mcparticle.pt(), mcparticle.phi() - phi);
        } else if (mcparticle.pdgCode() == 13) { // negative muon
          registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Neg"), mcparticle.pt(), mcparticle.phi() - phi);
        }
      }
    }
    return;
  }

  bool isSelectedMuon(const float pt, const float eta, const float rAtAbsorberEnd, const float pDCA, const float chi2, const uint8_t trackType, const float dcaXY)
  {
    if (pt < muoncuts.cfg_min_pt_track) {
      return false;
    }
    if (rAtAbsorberEnd < muoncuts.cfg_min_rabs_sa || muoncuts.cfg_max_rabs_sa < rAtAbsorberEnd) {
      return false;
    }
    if (rAtAbsorberEnd < muoncuts.cfg_mid_rabs ? pDCA > muoncuts.cfg_max_pdca_forSmallR : pDCA > muoncuts.cfg_max_pdca_forLargeR) {
      return false;
    }

    if (trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      if (eta < muoncuts.cfg_min_eta_track_gl || muoncuts.cfg_max_eta_track_gl < eta) {
        return false;
      }
      if (muoncuts.cfg_max_dcaxy_gl < dcaXY) {
        return false;
      }
      if (chi2 < 0.f || muoncuts.cfg_max_chi2_gl < chi2) {
        return false;
      }
      if (rAtAbsorberEnd < muoncuts.cfg_min_rabs_gl || muoncuts.cfg_max_rabs_gl < rAtAbsorberEnd) {
        return false;
      }
    } else if (trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      if (eta < muoncuts.cfg_min_eta_track_sa || muoncuts.cfg_max_eta_track_sa < eta) {
        return false;
      }
      if (chi2 < 0.f || muoncuts.cfg_max_chi2_sa < chi2) {
        return false;
      }
    } else {
      return false;
    }

    return true;
  }

  SliceCache cache;
  Preslice<aod::Tracks> perCollision_mid = o2::aod::track::collisionId;
  Preslice<aod::FwdTracks> perCollision_fwd = o2::aod::fwdtrack::collisionId;

  using MyCollisions = Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MyCollision = MyCollisions::iterator;

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
  using MyTrack = MyTracks::iterator;

  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  using MyFwdTrack = MyFwdTracks::iterator;

  using MyMFTTracks = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  using MyMFTTrack = MyMFTTracks::iterator;

  template <typename TCollision, typename TTrack>
  void fillElectron(TCollision const& collision, TTrack const& track, const float centrality)
  {
    if (cfgRequireGoodRCT && !rctCheckerCB.checkTable(collision)) {
      return;
    }
    auto mcparticle = track.template mcParticle_as<aod::McParticles>();

    if (std::abs(mcparticle.pdgCode()) != 11 || !(mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator())) {
      return;
    }
    if (cfg_reject_fake_match_its_tpc && o2::aod::pwgem::dilepton::utils::mcutil::hasFakeMatchITSTPC(track)) {
      return;
    }
    if (cfg_require_true_mc_collision_association && mcparticle.mcCollisionId() != collision.mcCollisionId()) {
      return;
    }

    if (!isSelectedTrack(track)) {
      return;
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    float pt = track_par_cov_recalc.getPt();
    float eta = track_par_cov_recalc.getEta();
    float phi = track_par_cov_recalc.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    if (!isSelectedTrackKine(track, pt, eta, dcaXY, dcaZ)) {
      return;
    }

    if (cfgFillTHnSparse) {
      registry.fill(HIST("Electron/hs_reso"), centrality, mcparticle.pt(), mcparticle.eta(), mcparticle.phi(), -mcparticle.pdgCode() / 11, (mcparticle.pt() - pt) / mcparticle.pt(), mcparticle.eta() - eta, mcparticle.phi() - phi);
    }
    if (cfgFillTH2) {
      registry.fill(HIST("Electron/hPt"), pt);
      registry.fill(HIST("Electron/hEtaPhi"), phi, eta);
      registry.fill(HIST("Electron/Ptgen_RelDeltaPt"), mcparticle.pt(), (mcparticle.pt() - pt) / mcparticle.pt());
      registry.fill(HIST("Electron/Ptgen_DeltaEta"), mcparticle.pt(), mcparticle.eta() - eta);
      if (mcparticle.pdgCode() == -11) { // positron
        registry.fill(HIST("Electron/Ptgen_DeltaPhi_Pos"), mcparticle.pt(), mcparticle.phi() - phi);
      } else if (mcparticle.pdgCode() == 11) { // electron
        registry.fill(HIST("Electron/Ptgen_DeltaPhi_Neg"), mcparticle.pt(), mcparticle.phi() - phi);
      }
    }
  }

  void processElectronSA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::McCollisions const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      // auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      // registry.fill(HIST("Event/Electron/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto tracks_per_coll = tracks.sliceBy(perCollision_mid, collision.globalIndex());
      for (const auto& track : tracks_per_coll) {
        if (!track.has_mcParticle()) {
          continue;
        }

        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        auto mccollision_from_mctrack = mctrack.template mcCollision_as<aod::McCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_mctrack.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }

        fillElectron(collision, track, centrality);
      } // end of track loop
    } // end of collision loop
  }
  PROCESS_SWITCH(CreateResolutionMap, processElectronSA, "create resolution map for electron at midrapidity", true);

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  void processElectronTTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const&, aod::TrackAssoc const& trackIndices, aod::McCollisions const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      // auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      // registry.fill(HIST("Event/Electron/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      for (const auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracks>();
        if (!track.has_mcParticle()) {
          continue;
        }
        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        auto mccollision_from_mctrack = mctrack.template mcCollision_as<aod::McCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_mctrack.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        fillElectron(collision, track, centrality);
      } // end of track loop
    } // end of collision loop
  }
  PROCESS_SWITCH(CreateResolutionMap, processElectronTTCA, "create resolution map for electron at midrapidity", false);

  Partition<MyFwdTracks> sa_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack); // MCH-MID
  Partition<MyFwdTracks> global_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack); // MFT-MCH-MID

  void processMuonSA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFwdTracks const&, MyMFTTracks const&, aod::McCollisions const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      // auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      // registry.fill(HIST("Event/Muon/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto sa_muons_per_coll = sa_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto global_muons_per_coll = global_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

      for (const auto& muon : sa_muons_per_coll) {
        if (!muon.has_mcParticle()) {
          continue;
        }

        auto mctrack = muon.template mcParticle_as<aod::McParticles>();
        auto mccollision_from_mctrack = mctrack.template mcCollision_as<aod::McCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_mctrack.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        fillMuon(collision, muon, centrality);
      } // end of standalone muon loop

      for (const auto& muon : global_muons_per_coll) {
        if (!muon.has_mcParticle()) {
          continue;
        }
        auto mctrack = muon.template mcParticle_as<aod::McParticles>();
        auto mccollision_from_mctrack = mctrack.template mcCollision_as<aod::McCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_mctrack.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        fillMuon(collision, muon, centrality);
      } // end of global muon loop

    } // end of collision loop
  }
  PROCESS_SWITCH(CreateResolutionMap, processMuonSA, "create resolution map for muon at forward rapidity", true);

  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  void processMuonTTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFwdTracks const&, MyMFTTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      // auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      // registry.fill(HIST("Event/Muon/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto muon = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (!muon.has_mcParticle()) {
          continue;
        }
        auto mctrack = muon.template mcParticle_as<aod::McParticles>();
        auto mccollision_from_mctrack = mctrack.template mcCollision_as<aod::McCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_mctrack.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        fillMuon(collision, muon, centrality);
      } // end of fwdtrack loop
    } // end of collision loop
  }
  PROCESS_SWITCH(CreateResolutionMap, processMuonTTCA, "create resolution map for muon at forward rapidity", false);

  void processGen(aod::McCollisions const& mcCollisions)
  {
    for (const auto& mccollision : mcCollisions) {
      registry.fill(HIST("Event/hGenID"), mccollision.getSubGeneratorId());
    }
  }
  PROCESS_SWITCH(CreateResolutionMap, processGen, "process generated info", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateResolutionMap>(cfgc, TaskName{"create-resolution-map"})};
}
