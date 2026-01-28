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
// Analysis task to check dca template for electrons/muons
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "TGeoGlobalMagField.h"

#include <array>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::fwdtrackutils;

struct checkMCTemplate {
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
  Configurable<uint> cfgDCAType{"cfgDCAType", 0, "type of DCA for output. 0:3D, 1:XY, 2:Z, else:3D"};

  ConfigurableAxis ConfPtBins{"ConfPtBins", {100, 0, 10}, "pT bins for output histograms"};
  ConfigurableAxis ConfEtaCBBins{"ConfEtaCBBins", {20, -1.0, +1.0}, "eta bins at midrapidity for output histograms"};
  ConfigurableAxis ConfEtaFWDBins{"ConfEtaFWDBins", {40, -4, -2}, "eta bins at forward rapidity for output histograms"};
  ConfigurableAxis ConfPhiBins{"ConfPhiBins", {36, 0, 2.f * M_PI}, "phi bins at forward rapidity for output histograms"};

  ConfigurableAxis ConfDCACBBins{"ConfDCACBBins", {1, 0, 10}, "dca bins for output histograms at midrapidity"};
  ConfigurableAxis ConfDCAFWDBins{"ConfDCAFWDBins", {1, 0, 10}, "dca bins for output histograms at fwd rapidity"};

  Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
  Configurable<std::string> cfgRCTLabelCB{"cfgRCTLabelCB", "CBT_hadronPID", "select 1 [CBT, CBT_hadron] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<std::string> cfgRCTLabelFWDGL{"cfgRCTLabelFWDGL", "CBT_muon_glo", "select 1 [CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
  Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
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
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
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
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<bool> cfg_reject_fake_match_its_tpc{"cfg_reject_fake_match_its_tpc", false, "flag to reject fake match between ITS-TPC"};
  } electroncuts;

  struct : ConfigurableGroup {
    std::string prefix = "muoncut_group";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
    Configurable<float> cfg_min_eta_track_gl{"cfg_min_eta_track_gl", -5.5, "min eta for global muon track"};
    Configurable<float> cfg_max_eta_track_gl{"cfg_max_eta_track_gl", -1.5, "max eta for global muon track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2_gl{"cfg_max_chi2_gl", 4, "max chi2/ndf for global muon track"};
    Configurable<float> cfg_max_chi2mft{"cfg_max_chi2mft", 1e+10, "max chi2/ndf for MFTsa track"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 40, "max chi2/ndf for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2/ndf for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy_gl{"cfg_max_dcaxy_gl", 0.1, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs_gl{"cfg_min_rabs_gl", 27.6, "min Radius at the absorber end for global muon track"};
    Configurable<float> cfg_max_rabs_gl{"cfg_max_rabs_gl", 89.5, "max Radius at the absorber end for global muon track"};
    Configurable<float> cfg_mid_rabs{"cfg_mid_rabs", 26.5, "middle R at absorber end for pDCA cut"};
    Configurable<float> cfg_max_pdca_forLargeR{"cfg_max_pdca_forLargeR", 324.f, "max. pDCA for large R at absorber end"};
    Configurable<float> cfg_max_pdca_forSmallR{"cfg_max_pdca_forSmallR", 594.f, "max. pDCA for small R at absorber end"};
    Configurable<float> cfg_max_reldpt{"cfg_max_reldpt", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_deta{"cfg_max_deta", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_dphi{"cfg_max_dphi", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
    Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
    Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "flag to require MFT hit map"};
    Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{4}, "hit map on MFT disks [0,1,2,3,4]. logical-OR of each double-sided disk"};
    Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
    Configurable<bool> cfg_reject_fake_match_mft_mch{"cfg_reject_fake_match_mft_mch", false, "flag to reject fake match between MFT-MCH"};
  } muoncuts;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::aod::rctsel::RCTFlagsChecker rctCheckerCB;
  o2::aod::rctsel::RCTFlagsChecker rctCheckerFWDGL;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBzMFT = 0;
  float d_bz = 0;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  ~checkMCTemplate() {}

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
    rctCheckerFWDGL.init(cfgRCTLabelFWDGL.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);

    mRunNumber = 0;
    d_bz = 0;
    mBzMFT = 0;

    const AxisSpec axis_pt{ConfPtBins, "p_{T,l} (GeV/c)"};
    const AxisSpec axis_eta_cb{ConfEtaCBBins, "#eta_{l}"};
    const AxisSpec axis_eta_fwd{ConfEtaFWDBins, "#eta_{l}"};
    const AxisSpec axis_phi{ConfPhiBins, "#varphi_{l} (rad.)"};
    const AxisSpec axis_dca_fwd{ConfDCAFWDBins, "DCA_{l}^{XY} (#sigma)"};
    const AxisSpec axis_charge{3, -1.5, +1.5, "sign"};

    std::string pair_dca_axis_title = "DCA_{l}^{3D} (#sigma)";
    if (cfgDCAType == 0) {
      pair_dca_axis_title = "DCA_{l}^{3D} (#sigma)";
    } else if (cfgDCAType == 1) {
      pair_dca_axis_title = "DCA_{l}^{XY} (#sigma)";
    } else if (cfgDCAType == 2) {
      pair_dca_axis_title = "DCA_{l}^{Z} (#sigma)";
    }
    const AxisSpec axis_dca_cb{ConfDCACBBins, pair_dca_axis_title};

    registry.add("Electron/c2l/hs", "hs", kTHnSparseF, {axis_pt, axis_eta_cb, axis_phi, axis_dca_cb, axis_charge}, false);
    registry.addClone("Electron/c2l/", "Electron/b2l/");
    registry.addClone("Electron/c2l/", "Electron/b2c2l/");

    registry.add("Muon/c2l/hs", "hs", kTHnSparseF, {axis_pt, axis_eta_fwd, axis_phi, axis_dca_fwd, axis_charge}, false);
    registry.addClone("Muon/c2l/", "Muon/b2l/");
    registry.addClone("Muon/c2l/", "Muon/b2c2l/");
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
      mBzMFT = d_bz;
      LOGF(info, "Bz at center of MFT = %f kZG manually", mBzMFT);
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

    // std::map<string, string> metadata;
    // auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    // auto ts = soreor.first;
    // auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    // o2::base::Propagator::initFieldFromGRP(grpmag);

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    const double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    mBzMFT = field->getBz(centerMFT); // Get field at centre of MFT
    LOGF(info, "Bz at center of MFT = %f kZG", mBzMFT);
    mRunNumber = bc.runNumber();
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
  bool isSelectedTrack(TTrack const& track, const float pt, const float eta, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
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

    if (track.tpcChi2NCl() > electroncuts.cfg_max_chi2tpc) {
      return false;
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

  template <bool withMFTCov, typename TCollision, typename TMuon, typename TMCParticles>
  void fillMuon(TCollision const& collision, TMuon const& muon, TMCParticles const& mcParticles)
  {
    auto mcparticle = muon.template mcParticle_as<aod::McParticles>();
    if (std::abs(mcparticle.pdgCode()) != 13 || !(mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) || !mcparticle.has_mothers()) {
      return;
    }
    if (cfg_require_true_mc_collision_association && mcparticle.mcCollisionId() != collision.mcCollisionId()) {
      return;
    }
    if (muon.trackType() != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      return;
    }

    if (std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(muon.globalIndex(), muon.matchMCHTrackId(), muon.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
      return;
    }

    if (muon.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return;
    }
    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(muon, muon, collision, propagationPoint::kToVertex, muoncuts.matchingZ, mBzMFT);

    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    float dcaX = propmuonAtPV.getX() - collision.posX();
    float dcaY = propmuonAtPV.getY() - collision.posY();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    float cXX = propmuonAtPV.getSigma2X();
    float cYY = propmuonAtPV.getSigma2Y();
    float cXY = propmuonAtPV.getSigmaXY();

    float dca = 999.f;
    float det = cXX * cYY - cXY * cXY; // determinanat
    if (det < 0) {
      dca = 999.f;
    } else {
      dca = std::sqrt(std::fabs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2.f * dcaX * dcaY * cXY) / det / 2.f)); // dca xy in sigma
    }

    // float rAtAbsorberEnd = muon.rAtAbsorberEnd(); // this works only for GlobalMuonTrack
    // float pDCA = propmuonAtPV.getP() * dcaXY;
    int nClustersMFT = 0;
    float ptMatchedMCHMID = propmuonAtPV.getPt();
    float etaMatchedMCHMID = propmuonAtPV.getEta();
    float phiMatchedMCHMID = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);

    // mcparticle for global MFT-MCH-MID is identical to mcparticle of MCH-MID track. If not, mismatch.
    auto mchtrack = muon.template matchMCHTrack_as<MyFwdTracks>(); // MCH-MID
    auto mfttrack = muon.template matchMFTTrack_as<MyMFTTracks>(); // MFTsa
    if (!mchtrack.has_mcParticle() || !mfttrack.has_mcParticle()) {
      return;
    }
    auto mcparticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>();
    auto mcparticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
    if (mcparticle.globalIndex() != mcparticle_MCHMID.globalIndex()) { // this should not happen. this is only for protection.
      return;
    }
    if (muoncuts.cfg_reject_fake_match_mft_mch && mcparticle.globalIndex() != mcparticle_MFT.globalIndex()) { // evaluate mismatch
      return;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToVertex, muoncuts.matchingZ, mBzMFT);
    ptMatchedMCHMID = propmuonAtPV_Matched.getPt();
    etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
    if (muoncuts.refitGlobalMuon) {
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
    }

    // o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToDCA, muoncuts.matchingZ, mBzMFT);
    // float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
    // float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
    // float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
    // pDCA = mchtrack.p() * dcaXY_Matched;
    // pDCA = propmuonAtPV.getP() * dcaXY;

    nClustersMFT = mfttrack.nClusters();
    float chi2mft = mfttrack.chi2() / (2.f * nClustersMFT - 5.f);

    if (nClustersMFT < muoncuts.cfg_min_ncluster_mft) {
      return;
    }

    if (chi2mft < 0.f || muoncuts.cfg_max_chi2mft < chi2mft) {
      return;
    }

    if (muon.chi2MatchMCHMFT() > muoncuts.cfg_max_matching_chi2_mftmch) {
      return;
    }

    float dpt = (ptMatchedMCHMID - pt) / pt;
    if (std::fabs(dpt) > muoncuts.cfg_max_reldpt) {
      return;
    }
    float deta = etaMatchedMCHMID - eta;
    float dphi = phiMatchedMCHMID - phi;
    o2::math_utils::bringToPMPi(dphi);
    if (std::sqrt(std::pow(deta / muoncuts.cfg_max_deta, 2) + std::pow(dphi / muoncuts.cfg_max_dphi, 2)) > 1.f) {
      return;
    }

    if (muoncuts.requireMFTHitMap) {
      std::vector<bool> hasMFTs{hasMFT<0, 1>(mfttrack), hasMFT<2, 3>(mfttrack), hasMFT<4, 5>(mfttrack), hasMFT<6, 7>(mfttrack), hasMFT<8, 9>(mfttrack)};
      for (int i = 0; i < static_cast<int>(muoncuts.requiredMFTDisks->size()); i++) {
        if (!hasMFTs[muoncuts.requiredMFTDisks->at(i)]) {
          return;
        }
      }
    }

    // fill histograms here
    auto mcmother = mcparticle.template mothers_as<aod::McParticles>()[0];

    if (isWeakDecayFromBeautyHadron(mcparticle, mcParticles)) { // hb->l is found in full decay chain.
      registry.fill(HIST("Muon/b2l/hs"), pt, eta, phi, dca, muon.sign());
    } else if (isWeakDecayFromCharmHadron(mcparticle, mcParticles)) { // hc->l is found in full decay chain.
      if (IsFromBeauty(mcmother, mcParticles) > 0) {                  // hb->hc->l is fond.
        registry.fill(HIST("Muon/b2c2l/hs"), pt, eta, phi, dca, muon.sign());
      } else { // prompt hc->l is found.
        registry.fill(HIST("Muon/c2l/hs"), pt, eta, phi, dca, muon.sign());
      }
    }
  }

  template <int begin = 0, int end = 9, typename T>
  bool hasMFT(T const& track)
  {
    // logical-OR
    uint64_t mftClusterSizesAndTrackFlags = track.mftClusterSizesAndTrackFlags();
    uint16_t clmap = 0;
    for (unsigned int layer = begin; layer <= end; layer++) {
      if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3f) {
        clmap |= (1 << layer);
      }
    }
    return (clmap > 0);
  }

  SliceCache cache;
  Preslice<aod::TracksIU> perCollision_mid = o2::aod::track::collisionId;
  Preslice<aod::FwdTracks> perCollision_fwd = o2::aod::fwdtrack::collisionId;
  PresliceUnsorted<aod::FwdTracks> fwdtracksPerMCHTrack = aod::fwdtrack::matchMCHTrackId;

  using MyCollisions = Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MyCollision = MyCollisions::iterator;

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
  using MyTrack = MyTracks::iterator;

  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  using MyFwdTrack = MyFwdTracks::iterator;

  using MyMFTTracks = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  using MyMFTTrack = MyMFTTracks::iterator;

  template <typename TCollision, typename TTrack, typename TMCParticles>
  void fillElectron(TCollision const& collision, TTrack const& track, TMCParticles const& mcParticles)
  {
    if (cfgRequireGoodRCT && !rctCheckerCB.checkTable(collision)) {
      return;
    }
    auto mcparticle = track.template mcParticle_as<aod::McParticles>();

    if (std::abs(mcparticle.pdgCode()) != 11 || !(mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) || !mcparticle.has_mothers()) {
      return;
    }
    if (electroncuts.cfg_reject_fake_match_its_tpc && o2::aod::pwgem::dilepton::utils::mcutil::hasFakeMatchITSTPC(track)) {
      return;
    }
    if (cfg_require_true_mc_collision_association && mcparticle.mcCollisionId() != collision.mcCollisionId()) {
      return;
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    float pt = trackParCov.getPt();
    float eta = trackParCov.getEta();
    float phi = trackParCov.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    if (!isSelectedTrack(track, pt, eta, dcaXY, dcaZ)) {
      return;
    }

    float dca_sigma = 999.f;
    if (cfgDCAType == 0) {
      float det = trackParCov.getSigmaY2() * trackParCov.getSigmaZ2() - trackParCov.getSigmaZY() * trackParCov.getSigmaZY();
      if (det < 0) {
        dca_sigma = 999.f;
      } else {
        float chi2 = (dcaXY * dcaXY * trackParCov.getSigmaZ2() + dcaZ * dcaZ * trackParCov.getSigmaY2() - 2. * dcaXY * dcaZ * trackParCov.getSigmaZY()) / det;
        dca_sigma = std::sqrt(std::fabs(chi2) / 2.);
      }
    } else if (cfgDCAType == 1) {
      dca_sigma = dcaXY / std::sqrt(trackParCov.getSigmaY2());
    } else if (cfgDCAType == 2) {
      dca_sigma = dcaZ / std::sqrt(trackParCov.getSigmaZ2());
    }

    // fill histograms here
    auto mcmother = mcparticle.template mothers_as<aod::McParticles>()[0];

    if (isWeakDecayFromBeautyHadron(mcparticle, mcParticles)) { // hb->l is found in full decay chain.
      registry.fill(HIST("Electron/b2l/hs"), pt, eta, phi, dca_sigma, track.sign());
    } else if (isWeakDecayFromCharmHadron(mcparticle, mcParticles)) { // hc->l is found in full decay chain.
      if (IsFromBeauty(mcmother, mcParticles) > 0) {                  // hb->hc->l is fond.
        registry.fill(HIST("Electron/b2c2l/hs"), pt, eta, phi, dca_sigma, track.sign());
      } else { // prompt hc->l is found.
        registry.fill(HIST("Electron/c2l/hs"), pt, eta, phi, dca_sigma, track.sign());
      }
    }
  }

  std::vector<std::tuple<int, int, int>> vec_min_chi2MatchMCHMFT; // std::pair<globalIndex of global muon, globalIndex of matched MCH-MID, globalIndex of MFT> -> chi2MatchMCHMFT;
  template <typename TMuons>
  void findBestMatchPerMCHMID(TMuons const& muons)
  {
    vec_min_chi2MatchMCHMFT.reserve(muons.size());
    for (const auto& muon : muons) {
      if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        const auto& muons_per_MCHMID = muons.sliceBy(fwdtracksPerMCHTrack, muon.globalIndex());
        // LOGF(info, "stanadalone: muon.globalIndex() = %d, muon.chi2MatchMCHMFT() = %f", muon.globalIndex(), muon.chi2MatchMCHMFT());
        // LOGF(info, "muons_per_MCHMID.size() = %d", muons_per_MCHMID.size());

        float min_chi2MatchMCHMFT = 1e+10;
        std::tuple<int, int, int> tupleIds_at_min;
        for (const auto& muon_tmp : muons_per_MCHMID) {
          if (muon_tmp.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
            // LOGF(info, "muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, muon_tmp.chi2MatchMCHMFT() = %f", muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId(), muon_tmp.chi2MatchMCHMFT());
            if (0.f < muon_tmp.chi2MatchMCHMFT() && muon_tmp.chi2MatchMCHMFT() < min_chi2MatchMCHMFT) {
              min_chi2MatchMCHMFT = muon_tmp.chi2MatchMCHMFT();
              tupleIds_at_min = std::make_tuple(muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId());
            }
          }
        }
        vec_min_chi2MatchMCHMFT.emplace_back(tupleIds_at_min);
        // LOGF(info, "min: muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, muon_tmp.chi2MatchMCHMFT() = %f", std::get<0>(tupleIds_at_min), std::get<1>(tupleIds_at_min), std::get<2>(tupleIds_at_min), min_chi2MatchMCHMFT);
      }
    } // end of muon loop
  }

  void processElectronSA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::McCollisions const&, aod::McParticles const& mcParticles)
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
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }

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

        fillElectron(collision, track, mcParticles);
      } // end of track loop
    } // end of collision loop
  }
  PROCESS_SWITCH(checkMCTemplate, processElectronSA, "check mc template for electron at midrapidity", true);

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  void processElectronTTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const&, aod::TrackAssoc const& trackIndices, aod::McCollisions const&, aod::McParticles const& mcParticles)
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
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }

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
        fillElectron(collision, track, mcParticles);
      } // end of track loop
    } // end of collision loop
  }
  PROCESS_SWITCH(checkMCTemplate, processElectronTTCA, "check mc template for electron at midrapidity", false);

  void processMuonSA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFwdTracks const& fwdtracks, MyMFTTracks const&, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    findBestMatchPerMCHMID(fwdtracks);
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
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }

      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision_fwd, collision.globalIndex());
      for (const auto& muon : fwdtracks_per_coll) {
        if (!muon.has_mcParticle()) {
          continue;
        }

        auto mctrack = muon.template mcParticle_as<aod::McParticles>();
        auto mccollision_from_mctrack = mctrack.template mcCollision_as<aod::McCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_mctrack.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        fillMuon<false>(collision, muon, mcParticles);
      } // end of standalone muon loop

    } // end of collision loop
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
  }
  PROCESS_SWITCH(checkMCTemplate, processMuonSA, "check mc template for muon at forward rapidity", false);

  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  void processMuonTTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFwdTracks const& fwdtracks, MyMFTTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    findBestMatchPerMCHMID(fwdtracks);
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
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }

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
        fillMuon<false>(collision, muon, mcParticles);
      } // end of fwdtrack loop
    } // end of collision loop
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
  }
  PROCESS_SWITCH(checkMCTemplate, processMuonTTCA, "check mc template for muon at forward rapidity", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<checkMCTemplate>(cfgc, TaskName{"check-mc-template"})};
}
