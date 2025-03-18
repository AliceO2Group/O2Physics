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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

struct CreateResolutionMap {
  // Index used to set different options for Muon propagation
  enum class MuonExtrapolation : int {
    kToVertex = 0, // propagtion to vertex by default
    kToDCA = 1,
    kToRabs = 2,
  };
  using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<double, 5>;

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
  Configurable<bool> cfg_require_true_mc_collision_association{"cfg_require_true_mc_collision_association", true, "flag to require true mc collision association"};
  Configurable<bool> cfg_reject_fake_match_its_tpc{"cfg_reject_fake_match_its_tpc", false, "flag to reject fake match between ITS-TPC"};
  // Configurable<bool> cfg_reject_fake_match_its_tpc_tof{"cfg_reject_fake_match_its_tpc_tof", false, "flag to reject fake match between ITS-TPC-TOF"};
  Configurable<bool> cfg_reject_fake_match_mft_mch{"cfg_reject_fake_match_mft_mch", false, "flag to reject fake match between MFT-MCH"};

  ConfigurableAxis ConfPtGenBins{"ConfPtGenBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00}, "gen. pT bins for output histograms"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0, 10, 30, 50, 110}, "centrality (%) bins for output histograms"};

  ConfigurableAxis ConfEtaCBGenBins{"ConfEtaCBGenBins", {30, -1.5, +1.5}, "gen. eta bins at midrapidity for output histograms"};
  ConfigurableAxis ConfEtaFWDGenBins{"ConfEtaFWDGenBins", {40, -5.5, -1.5}, "gen. eta bins at forward rapidity for output histograms"};
  ConfigurableAxis ConfPhiGenBins{"ConfPhiGenBins", {72, 0, 2.f * M_PI}, "gen. eta bins at forward rapidity for output histograms"};

  ConfigurableAxis ConfRelDeltaPtBins{"ConfRelDeltaPtBins", {200, -1.f, +1.f}, "rel. dpt for output histograms"};
  ConfigurableAxis ConfDeltaEtaBins{"ConfDeltaEtaBins", {100, -0.1f, +0.1f}, "deta bins for output histograms"};
  ConfigurableAxis ConfDeltaPhiBins{"ConfDeltaPhiBins", {100, -0.1f, +0.1f}, "dphi bins for output histograms"};

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
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 90, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 0, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_tpc_cr_findable_ratio{"cfg_min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.2, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.2, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
  } electroncuts;

  struct : ConfigurableGroup {
    std::string prefix = "muoncut_group";
    Configurable<uint> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -5.5, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -1.5, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+10, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 1e+10, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
  } muoncuts;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::globaltracking::MatchGlobalFwd mMatching;
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

    registry.add("Event/Electron/hImpPar_Centrality", "true imapact parameter vs. estimated centrality;impact parameter (fm);centrality (%)", kTH2F, {{200, 0, 20}, {110, 0, 110}}, true);
    registry.add("Event/Muon/hImpPar_Centrality", "true imapact parameter vs. estimated centrality;impact parameter (fm);centrality (%)", kTH2F, {{200, 0, 20}, {110, 0, 110}}, true);

    registry.add("Electron/hPt", "rec. p_{T,l};p_{T,l} (GeV/c)", kTH1F, {{1000, 0, 10}}, false);
    registry.add("Electron/hEtaPhi", "rec. #eta vs. #varphi;#varphi_{l} (rad.);#eta_{l}", kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    registry.add("Electron/Ptgen_RelDeltaPt", "resolution", kTH2F, {{axis_pt_gen}, {axis_dpt}}, true);
    registry.add("Electron/Ptgen_DeltaEta", "resolution", kTH2F, {{axis_pt_gen}, {axis_deta}}, true);
    registry.add("Electron/Ptgen_DeltaPhi_Pos", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
    registry.add("Electron/Ptgen_DeltaPhi_Neg", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
    registry.addClone("Electron/", "StandaloneMuon/");
    registry.addClone("Electron/", "GlobalMuon/");

    registry.add("Electron/hs_reso", "8D resolution positive", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_cb_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
    registry.add("StandaloneMuon/hs_reso", "8D resolution positive", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
    registry.add("GlobalMuon/hs_reso", "8D resolution positive", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
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

  template <typename TCollision, typename TTrack>
  bool isSelectedTrack(TCollision const& collision, TTrack const& track)
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

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    // LOGF(info, "collision.globalIndex() = %d, track.collisionId() = %d, track.pt() = %.16f, track_par_cov_recalc.getPt() = %.16f", collision.globalIndex(), track.collisionId(), track.pt(), track_par_cov_recalc.getPt());

    if (std::fabs(dcaXY) > electroncuts.cfg_max_dcaxy || std::fabs(dcaZ) > electroncuts.cfg_max_dcaz) {
      return false;
    }

    if (track_par_cov_recalc.getPt() < electroncuts.cfg_min_pt_track || std::fabs(track_par_cov_recalc.getEta()) > electroncuts.cfg_max_eta_track) {
      return false;
    }

    return true;
  }

  template <typename T, typename C>
  o2::dataformats::GlobalFwdTrack PropagateMuon(T const& muon, C const& collision, const CreateResolutionMap::MuonExtrapolation endPoint)
  {
    double chi2 = muon.chi2();
    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<float> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                          muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                          muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack propmuon;

    if (static_cast<int>(muon.trackType()) > 2) { // MCH-MID or MCH standalone
      o2::dataformats::GlobalFwdTrack track;
      track.setParameters(tpars);
      track.setZ(fwdtrack.getZ());
      track.setCovariances(tcovs);
      auto mchTrack = mMatching.FwdtoMCH(track);

      if (endPoint == CreateResolutionMap::MuonExtrapolation::kToVertex) {
        o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
      }
      if (endPoint == CreateResolutionMap::MuonExtrapolation::kToDCA) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
      }
      if (endPoint == CreateResolutionMap::MuonExtrapolation::kToRabs) {
        o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
      }

      auto proptrack = mMatching.MCHtoFwd(mchTrack);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());
    } else if (static_cast<int>(muon.trackType()) < 2) { // MFT-MCH-MID
      double centerMFT[3] = {0, 0, -61.4};
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
      auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.x(), muon.y(), muon.z(), collision.posX(), collision.posY(), collision.posZ());
      auto x2x0 = static_cast<float>(geoMan.meanX2X0);
      fwdtrack.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x0);
      propmuon.setParameters(fwdtrack.getParameters());
      propmuon.setZ(fwdtrack.getZ());
      propmuon.setCovariances(fwdtrack.getCovariances());
    }

    v1.clear();
    v1.shrink_to_fit();

    return propmuon;
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

    o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMuon(muon, collision, CreateResolutionMap::MuonExtrapolation::kToVertex);
    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();

    if (pt < muoncuts.cfg_min_pt_track) {
      return;
    }

    if (eta < muoncuts.cfg_min_eta_track || muoncuts.cfg_max_eta_track < eta) {
      return;
    }

    o2::math_utils::bringTo02Pi(phi);
    if (phi < 0.f || 2.f * M_PI < phi) {
      return;
    }

    float rAtAbsorberEnd = muon.rAtAbsorberEnd();
    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(muon, collision, CreateResolutionMap::MuonExtrapolation::kToRabs);
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    }

    if (rAtAbsorberEnd < muoncuts.cfg_min_rabs || muoncuts.cfg_max_rabs < rAtAbsorberEnd) {
      return;
    }

    if (rAtAbsorberEnd < 26.5) {
      if (muon.pDca() > 594.f) {
        return;
      }
    } else {
      if (muon.pDca() > 324.f) {
        return;
      }
    }

    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && muon.chi2MatchMCHMFT() > muoncuts.cfg_max_matching_chi2_mftmch) {
      return;
    }

    if (cfg_reject_fake_match_mft_mch && muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && o2::aod::pwgem::dilepton::utils::mcutil::hasFakeMatchMFTMCH(muon)) {
      return;
    }

    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      registry.fill(HIST("StandaloneMuon/hPt"), pt);
      registry.fill(HIST("StandaloneMuon/hEtaPhi"), phi, eta);
      registry.fill(HIST("StandaloneMuon/hs_reso"), centrality, mcparticle.pt(), mcparticle.eta(), mcparticle.phi(), -mcparticle.pdgCode() / 13, (mcparticle.pt() - pt) / mcparticle.pt(), mcparticle.eta() - eta, mcparticle.phi() - phi);
      registry.fill(HIST("StandaloneMuon/Ptgen_RelDeltaPt"), mcparticle.pt(), (mcparticle.pt() - pt) / mcparticle.pt());
      registry.fill(HIST("StandaloneMuon/Ptgen_DeltaEta"), mcparticle.pt(), mcparticle.eta() - eta);
      if (mcparticle.pdgCode() == -13) { // positive muon
        registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Pos"), mcparticle.pt(), mcparticle.phi() - phi);
      } else if (mcparticle.pdgCode() == 13) { // negative muon
        registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Neg"), mcparticle.pt(), mcparticle.phi() - phi);
      }
    } else if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      registry.fill(HIST("GlobalMuon/hPt"), pt);
      registry.fill(HIST("GlobalMuon/hEtaPhi"), phi, eta);
      registry.fill(HIST("GlobalMuon/hs_reso"), centrality, mcparticle.pt(), mcparticle.eta(), mcparticle.phi(), -mcparticle.pdgCode() / 13, (mcparticle.pt() - pt) / mcparticle.pt(), mcparticle.eta() - eta, mcparticle.phi() - phi);
      registry.fill(HIST("GlobalMuon/Ptgen_RelDeltaPt"), mcparticle.pt(), (mcparticle.pt() - pt) / mcparticle.pt());
      registry.fill(HIST("GlobalMuon/Ptgen_DeltaEta"), mcparticle.pt(), mcparticle.eta() - eta);
      if (mcparticle.pdgCode() == -13) { // positive muon
        registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Pos"), mcparticle.pt(), mcparticle.phi() - phi);
      } else if (mcparticle.pdgCode() == 13) { // negative muon
        registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Neg"), mcparticle.pt(), mcparticle.phi() - phi);
      }
    }
    return;
  }

  SliceCache cache;
  Preslice<aod::Tracks> perCollision_mid = o2::aod::track::collisionId;
  Preslice<aod::FwdTracks> perCollision_fwd = o2::aod::fwdtrack::collisionId;

  using MyCollisions = Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MyCollision = MyCollisions::iterator;

  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov, aod::McTrackLabels>;
  using MyTrack = MyTracks::iterator;

  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  using MyFwdTrack = MyFwdTracks::iterator;

  template <typename TCollision, typename TTrack>
  void fillElectron(TCollision const& collision, TTrack const& track, const float centrality)
  {
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
    if (!isSelectedTrack(collision, track)) {
      return;
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc, 2.f, matCorr, &mDcaInfoCov);
    // float dcaXY = mDcaInfoCov.getY();
    // float dcaZ = mDcaInfoCov.getZ();

    float pt = track_par_cov_recalc.getPt();
    float eta = track_par_cov_recalc.getEta();
    float phi = track_par_cov_recalc.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    registry.fill(HIST("Electron/hPt"), pt);
    registry.fill(HIST("Electron/hEtaPhi"), phi, eta);
    registry.fill(HIST("Electron/hs_reso"), centrality, mcparticle.pt(), mcparticle.eta(), mcparticle.phi(), -mcparticle.pdgCode() / 11, (mcparticle.pt() - pt) / mcparticle.pt(), mcparticle.eta() - eta, mcparticle.phi() - phi);
    registry.fill(HIST("Electron/Ptgen_RelDeltaPt"), mcparticle.pt(), (mcparticle.pt() - pt) / mcparticle.pt());
    registry.fill(HIST("Electron/Ptgen_DeltaEta"), mcparticle.pt(), mcparticle.eta() - eta);
    if (mcparticle.pdgCode() == -11) { // positron
      registry.fill(HIST("Electron/Ptgen_DeltaPhi_Pos"), mcparticle.pt(), mcparticle.phi() - phi);
    } else if (mcparticle.pdgCode() == 11) { // electron
      registry.fill(HIST("Electron/Ptgen_DeltaPhi_Neg"), mcparticle.pt(), mcparticle.phi() - phi);
    }
  }

  void processElectronSA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::McCollisions const&, aod::McParticles const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      registry.fill(HIST("Event/Electron/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto tracks_per_coll = tracks.sliceBy(perCollision_mid, collision.globalIndex());
      for (auto& track : tracks_per_coll) {
        if (!track.has_mcParticle()) {
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
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      registry.fill(HIST("Event/Electron/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracks>();
        if (!track.has_mcParticle()) {
          continue;
        }
        fillElectron(collision, track, centrality);
      } // end of track loop
    } // end of collision loop
  }
  PROCESS_SWITCH(CreateResolutionMap, processElectronTTCA, "create resolution map for electron at midrapidity", false);

  Partition<MyFwdTracks> sa_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack); // MCH-MID
  Partition<MyFwdTracks> global_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack); // MFT-MCH-MID

  void processMuonSA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFwdTracks const&, aod::McCollisions const&, aod::McParticles const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      registry.fill(HIST("Event/Muon/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto sa_muons_per_coll = sa_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto global_muons_per_coll = global_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

      for (auto& muon : sa_muons_per_coll) {
        if (!muon.has_mcParticle()) {
          continue;
        }
        fillMuon(collision, muon, centrality);
      } // end of standalone muon loop

      for (auto& muon : global_muons_per_coll) {
        if (!muon.has_mcParticle()) {
          continue;
        }
        fillMuon(collision, muon, centrality);
      } // end of global muon loop

    } // end of collision loop
  }
  PROCESS_SWITCH(CreateResolutionMap, processMuonSA, "create resolution map for muon at forward rapidity", true);

  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  void processMuonTTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFwdTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!isSelectedEvent(collision)) {
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      registry.fill(HIST("Event/Muon/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto muon = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (!muon.has_mcParticle()) {
          continue;
        }
        fillMuon(collision, muon, centrality);
      } // end of fwdtrack loop
    } // end of collision loop
  }
  PROCESS_SWITCH(CreateResolutionMap, processMuonTTCA, "create resolution map for muon at forward rapidity", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateResolutionMap>(cfgc, TaskName{"create-resolution-map"})};
}
