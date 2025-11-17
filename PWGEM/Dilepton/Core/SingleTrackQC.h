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
// This code runs loop over electrons for QC.
//    Please write to: daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_CORE_SINGLETRACKQC_H_
#define PWGEM_DILEPTON_CORE_SINGLETRACKQC_H_

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/MlResponseDielectronSingleTrack.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Tools/ML/MlResponse.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TString.h"

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBits>;
using MyCollisionWithSWT = MyCollisionsWithSWT::iterator;

using MyElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyElectron = MyElectrons::iterator;
using FilteredMyElectrons = soa::Filtered<MyElectrons>;
using FilteredMyElectron = FilteredMyElectrons::iterator;

using MyMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMGlobalMuonSelfIds>;
using MyMuon = MyMuons::iterator;
using FilteredMyMuons = soa::Filtered<MyMuons>;
using FilteredMyMuon = FilteredMyMuons::iterator;

template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename... Types>
struct SingleTrackQC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<std::string> cfg_swt_name{"cfg_swt_name", "fHighTrackMult", "desired software trigger name"}; // 1 trigger name per 1 task. fHighTrackMult, fHighFt0Mult
  Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
  Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
  Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};

  ConfigurableAxis ConfPtlBins{"ConfPtlBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTl bins for output histograms"};
  ConfigurableAxis ConfDCA3DBins{"ConfDCA3DBins", {VARIABLE_WIDTH, 0.0, 10.0}, "DCA3d bins in sigma for output histograms"};
  ConfigurableAxis ConfDCAXYBins{"ConfDCAXYBins", {VARIABLE_WIDTH, -10.0, 10.0}, "DCAxy bins in sigma for output histograms"};
  ConfigurableAxis ConfDCAZBins{"ConfDCAZBins", {VARIABLE_WIDTH, -10.0, 10.0}, "DCAz bins in sigma for output histograms"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
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
    Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "number of inactive chips on ITS layer 3 are below threshold "};
    Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "number of inactive chips on ITS layers 0-3 are below threshold "};
    Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "number of inactive chips on all ITS layers are below threshold "};
    // for RCT
    Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } eventcuts;

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "dielectroncut_group";

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<bool> cfg_mirror_phi_track{"cfg_mirror_phi_track", false, "mirror the phi cut around Pi, min and max phi should be in 0-Pi"};
    Configurable<bool> cfg_reject_phi_track{"cfg_reject_phi_track", false, "reject the phi interval"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2 TOF"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.2, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.2, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_min_its_cluster_size{"cfg_min_its_cluster_size", 0.f, "min ITS cluster size"};
    Configurable<float> cfg_max_its_cluster_size{"cfg_max_its_cluster_size", 16.f, "max ITS cluster size"};
    Configurable<float> cfg_min_rel_diff_pin{"cfg_min_rel_diff_pin", -1e+10, "min rel. diff. between pin and ppv"};
    Configurable<float> cfg_max_rel_diff_pin{"cfg_max_rel_diff_pin", +1e+10, "max rel. diff. between pin and ppv"};
    Configurable<float> cfgRefR{"cfgRefR", 0.50, "ref. radius (m) for calculating phi position"}; // 0.50 +/- 0.06 can be syst. unc.
    Configurable<float> cfg_min_phiposition_track{"cfg_min_phiposition_track", 0.f, "min phi position for single track at certain radius"};
    Configurable<float> cfg_max_phiposition_track{"cfg_max_phiposition_track", 6.3, "max phi position for single track at certain radius"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3, kTOFif = 4, kPIDML = 5]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    // Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    // Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_min_pin_pirejTPC{"cfg_min_pin_pirejTPC", 0.f, "min. pin for pion rejection in TPC"};
    Configurable<float> cfg_max_pin_pirejTPC{"cfg_max_pin_pirejTPC", 1e+10, "max. pin for pion rejection in TPC"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
    Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to enable ITSsa tracks"};
    Configurable<float> cfg_max_pt_track_ITSsa{"cfg_max_pt_track_ITSsa", 0.15, "max pt for ITSsa tracks"};

    // configuration for PID ML
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
    Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
    Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{-999999., 999999.}, "Bin limits for ML application"};
    Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.95}, "ML cuts per bin"};
    Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature"}, "Names of ML model input features"};
    Configurable<std::string> nameBinningFeature{"nameBinningFeature", "pt", "Names of ML model binning feature"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dielectroncuts;

  DimuonCut fDimuonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dimuoncut_group";

    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "max phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+6, "max chi2/ndf"};
    Configurable<float> cfg_max_chi2mft{"cfg_max_chi2mft", 1e+6, "max chi2/ndf"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 40, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
    Configurable<float> cfg_max_relDPt_wrt_matchedMCHMID{"cfg_max_relDPt_wrt_matchedMCHMID", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_DEta_wrt_matchedMCHMID{"cfg_max_DEta_wrt_matchedMCHMID", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_DPhi_wrt_matchedMCHMID{"cfg_max_DPhi_wrt_matchedMCHMID", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
    Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "flag to apply MFT hit map"};
    Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{0}, "hit map on MFT disks [0,1,2,3,4]. logical-OR of each double-sided disk"};
  } dimuoncuts;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false}; // 1 HistogramRegistry can keep up to 512 histograms
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};

  ~SingleTrackQC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      const AxisSpec axis_pt{ConfPtlBins, "p_{T,e} (GeV/c)"};
      const AxisSpec axis_eta{20, -1.0, +1.0, "#eta_{e}"};
      const AxisSpec axis_phi{36, 0.0, 2 * M_PI, "#varphi_{e} (rad.)"};
      const AxisSpec axis_phiposition{36, 0.0, 2 * M_PI, "#varphi_{e}^{*} (rad.)"};
      const AxisSpec axis_dca3D{ConfDCA3DBins, "DCA_{e}^{3D} (#sigma)"};
      const AxisSpec axis_dcaXY{ConfDCAXYBins, "DCA_{e}^{XY} (#sigma)"};
      const AxisSpec axis_dcaZ{ConfDCAZBins, "DCA_{e}^{Z} (#sigma)"};

      // track info
      fRegistry.add("Track/positive/hs", "rec. single electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_dca3D, axis_dcaXY, axis_dcaZ}, true);
      fRegistry.add("Track/positive/hPhiPosition", Form("phi position at r_{xy} = %3.2f m", dielectroncuts.cfgRefR.value), kTH1F, {axis_phiposition}, false);
      fRegistry.add("Track/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{4000, -20, 20}}, false);
      fRegistry.add("Track/positive/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.f, 1.f}}, false);
      fRegistry.add("Track/positive/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{400, -20.0f, 20.0f}, {400, -20.0f, 20.0f}}, false);
      fRegistry.add("Track/positive/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/positive/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/positive/hDCA3dRes_Pt", "DCA_{3D} resolution vs. pT;p_{T} (GeV/c);DCA_{3D} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/positive/hNclsTPC_Pt", "number of TPC clusters;p_{T,e} (GeV/c);TPC N_{cls}", kTH2F, {axis_pt, {161, -0.5, 160.5}}, false);
      fRegistry.add("Track/positive/hNcrTPC_Pt", "number of TPC crossed rows;p_{T,e} (GeV/c);TPC N_{CR}", kTH2F, {axis_pt, {161, -0.5, 160.5}}, false);
      fRegistry.add("Track/positive/hChi2TPC", "chi2/number of TPC clusters;TPC #chi^{2}/N_{CR}", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/positive/hDeltaPin", "p_{in} vs. p_{pv};p_{in} (GeV/c);(p_{pv} - p_{in})/p_{in}", kTH2F, {{1000, 0, 10}, {200, -1, +1}}, false);
      fRegistry.add("Track/positive/hTPCNcr2Nf", "TPC Ncr/Nfindable;TPC N_{CR}/N_{cls}^{findable}", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/positive/hTPCNcls2Nf", "TPC Ncls/Nfindable;TPC N_{cls}/N_{cls}^{findable}", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/positive/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("Track/positive/hNclsITS", "number of ITS clusters;ITS N_{cls}", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("Track/positive/hChi2ITS", "chi2/number of ITS clusters;ITS #chi^{2}/N_{cls}", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/positive/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("Track/positive/hChi2TOF", "TOF Chi2;p_{pv} (GeV/c);TOF #chi^{2}", kTH2F, {{1000, 0, 10}, {100, 0, 10}}, false);

      fRegistry.add("Track/positive/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      // fRegistry.add("Track/positive/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

      fRegistry.add("Track/positive/hTOFbeta", "TOF #beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
      fRegistry.add("Track/positive/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      // fRegistry.add("Track/positive/hTOFNsigmaMu", "TOF n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      // fRegistry.add("Track/positive/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      // fRegistry.add("Track/positive/hTOFNsigmaKa", "TOF n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      // fRegistry.add("Track/positive/hTOFNsigmaPr", "TOF n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

      fRegistry.add("Track/positive/hProbElBDT", "probability to be e from BDT;p_{in} (GeV/c);BDT score;", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("Track/positive/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda);", kTH2F, {{1000, 0.f, 10.f}, {150, 0, 15}}, false);
      fRegistry.add("Track/positive/hMeanClusterSizeITSib", "mean cluster size ITS inner barrel;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda);", kTH2F, {{1000, 0.f, 10.f}, {150, 0, 15}}, false);
      fRegistry.add("Track/positive/hMeanClusterSizeITSob", "mean cluster size ITS outer barrel;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda);", kTH2F, {{1000, 0.f, 10.f}, {150, 0, 15}}, false);

      fRegistry.addClone("Track/positive/", "Track/negative/");
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      const AxisSpec axis_pt{ConfPtlBins, "p_{T,#mu} (GeV/c)"};
      const AxisSpec axis_eta{50, -6, -1, "#eta_{#mu}"};
      const AxisSpec axis_phi{36, 0, 2 * M_PI, "#varphi_{#mu} (rad.)"};
      const AxisSpec axis_dca{ConfDCAXYBins, "DCA_{#mu}^{XY} (#sigma)"};

      // track info
      fRegistry.add("Track/positive/hs", "rec. single muon", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_dca}, true);
      fRegistry.add("Track/positive/hEtaPhi_MatchMCHMID", "#eta vs. #varphi of matched MCHMID", kTH2F, {{180, 0, 2.f * M_PI}, {100, -6, -1}}, false);
      fRegistry.add("Track/positive/hsDelta", "diff. between GL and associated SA;p_{T}^{gl} (GeV/c);(p_{T}^{sa} - p_{T}^{gl})/p_{T}^{gl};#Delta#eta;#Delta#varphi (rad.);", kTHnSparseF, {axis_pt, {100, -0.5, +0.5}, {100, -0.5, +0.5}, {90, -M_PI / 4, M_PI / 4}}, false);
      fRegistry.add("Track/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -5, 5}}, false);
      fRegistry.add("Track/positive/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
      fRegistry.add("Track/positive/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -0.5f, 0.5f}, {200, -0.5f, 0.5f}}, false);
      fRegistry.add("Track/positive/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/positive/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0, 500}}, false);
      fRegistry.add("Track/positive/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0, 500}}, false);
      fRegistry.add("Track/positive/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0, 500}}, false);
      fRegistry.add("Track/positive/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
      fRegistry.add("Track/positive/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
      fRegistry.add("Track/positive/hPDCA", "pDCA;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 100}, {100, 0.0f, 1000}}, false);
      fRegistry.add("Track/positive/hChi2", "chi2;chi2/ndf", kTH1F, {{100, 0.0f, 10}}, false);
      fRegistry.add("Track/positive/hChi2MFT", "chi2MFT;chi2/ndf", kTH1F, {{100, 0.0f, 10}}, false);
      fRegistry.add("Track/positive/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/positive/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/positive/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
      fRegistry.addClone("Track/positive/", "Track/negative/");
    }
  }

  int mRunNumber;
  float d_bz;
  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);

    DefineEMEventCut();
    DefineDielectronCut();
    DefineDimuonCut();
    addhistograms();
    mRunNumber = 0;
    d_bz = 0;

    if (doprocessNorm) {
      fRegistry.addClone("Event/before/hCollisionCounter", "Event/norm/hCollisionCounter");
    }
    if (doprocessQC_TriggeredData) {
      LOGF(info, "Trigger analysis is enabled. Desired trigger name = %s", cfg_swt_name.value.data());
      fRegistry.add("NormTrigger/hInspectedTVX", "inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
      fRegistry.add("NormTrigger/hScalers", "trigger counter before DS;run number;counter", kTProfile, {{80000, 520000.5, 600000.5}}, true);
      fRegistry.add("NormTrigger/hSelections", "trigger counter after DS;run number;counter", kTProfile, {{80000, 520000.5, 600000.5}}, true);
      auto hTriggerCounter = fRegistry.add<TH2>("NormTrigger/hTriggerCounter", Form("trigger counter of %s;run number;", cfg_swt_name.value.data()), kTH2D, {{80000, 520000.5, 600000.5}, {2, -0.5, 1.5}}, false);
      hTriggerCounter->GetYaxis()->SetBinLabel(1, "Analyzed Trigger");
      hTriggerCounter->GetYaxis()->SetBinLabel(2, "Analyzed TOI");
    }
    if (doprocessBC) {
      auto hTVXCounter = fRegistry.add<TH1>("BC/hTVXCounter", "TVX counter", kTH1D, {{6, -0.5f, 5.5f}});
      hTVXCounter->GetXaxis()->SetBinLabel(1, "TVX");
      hTVXCounter->GetXaxis()->SetBinLabel(2, "TVX && NoTFB");
      hTVXCounter->GetXaxis()->SetBinLabel(3, "TVX && NoITSROFB");
      hTVXCounter->GetXaxis()->SetBinLabel(4, "TVX && GoodRCT");
      hTVXCounter->GetXaxis()->SetBinLabel(5, "TVX && NoTFB && NoITSROFB");
      hTVXCounter->GetXaxis()->SetBinLabel(6, "TVX && NoTFB && NoITSROFB && GoodRCT");
    }
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kG";
    }

    mRunNumber = collision.runNumber();
    fDielectronCut.SetTrackPhiPositionRange(dielectroncuts.cfg_min_phiposition_track, dielectroncuts.cfg_max_phiposition_track, dielectroncuts.cfgRefR, d_bz, dielectroncuts.cfg_mirror_phi_track);
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(eventcuts.cfgZvtxMin, eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireVertexTOFmatched(eventcuts.cfgRequireVertexTOFmatched);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
    fEMEventCut.SetRequireGoodITSLayer3(eventcuts.cfgRequireGoodITSLayer3);
    fEMEventCut.SetRequireGoodITSLayer0123(eventcuts.cfgRequireGoodITSLayer0123);
    fEMEventCut.SetRequireGoodITSLayersAll(eventcuts.cfgRequireGoodITSLayersAll);
  }

  o2::analysis::MlResponseDielectronSingleTrack<float> mlResponseSingleTrack;
  void DefineDielectronCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, dielectroncuts.cfg_max_pt_track);
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetTrackPhiRange(dielectroncuts.cfg_min_phi_track, dielectroncuts.cfg_max_phi_track, dielectroncuts.cfg_mirror_phi_track, dielectroncuts.cfg_reject_phi_track);
    fDielectronCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetMaxFracSharedClustersTPC(dielectroncuts.cfg_max_frac_shared_clusters_tpc);
    fDielectronCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITS(dielectroncuts.cfg_min_its_cluster_size, dielectroncuts.cfg_max_its_cluster_size);
    fDielectronCut.SetTrackMaxDcaXY(dielectroncuts.cfg_max_dcaxy);
    fDielectronCut.SetTrackMaxDcaZ(dielectroncuts.cfg_max_dcaz);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);
    fDielectronCut.SetChi2TOF(0.0, dielectroncuts.cfg_max_chi2tof);
    fDielectronCut.SetRelDiffPin(dielectroncuts.cfg_min_rel_diff_pin, dielectroncuts.cfg_max_rel_diff_pin);
    fDielectronCut.IncludeITSsa(dielectroncuts.includeITSsa, dielectroncuts.cfg_max_pt_track_ITSsa);

    // for eID
    fDielectronCut.SetPIDScheme(dielectroncuts.cfg_pid_scheme);
    fDielectronCut.SetTPCNsigmaElRange(dielectroncuts.cfg_min_TPCNsigmaEl, dielectroncuts.cfg_max_TPCNsigmaEl);
    // fDielectronCut.SetTPCNsigmaMuRange(dielectroncuts.cfg_min_TPCNsigmaMu, dielectroncuts.cfg_max_TPCNsigmaMu);
    fDielectronCut.SetTPCNsigmaPiRange(dielectroncuts.cfg_min_TPCNsigmaPi, dielectroncuts.cfg_max_TPCNsigmaPi);
    fDielectronCut.SetTPCNsigmaKaRange(dielectroncuts.cfg_min_TPCNsigmaKa, dielectroncuts.cfg_max_TPCNsigmaKa);
    fDielectronCut.SetTPCNsigmaPrRange(dielectroncuts.cfg_min_TPCNsigmaPr, dielectroncuts.cfg_max_TPCNsigmaPr);
    fDielectronCut.SetTOFNsigmaElRange(dielectroncuts.cfg_min_TOFNsigmaEl, dielectroncuts.cfg_max_TOFNsigmaEl);
    fDielectronCut.SetPinRangeForPionRejectionTPC(dielectroncuts.cfg_min_pin_pirejTPC, dielectroncuts.cfg_max_pin_pirejTPC);

    if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      std::vector<float> binsML{};
      binsML.reserve(dielectroncuts.binsMl.value.size());
      for (size_t i = 0; i < dielectroncuts.binsMl.value.size(); i++) {
        binsML.emplace_back(dielectroncuts.binsMl.value[i]);
      }
      std::vector<float> thresholdsML{};
      thresholdsML.reserve(dielectroncuts.cutsMl.value.size());
      for (size_t i = 0; i < dielectroncuts.cutsMl.value.size(); i++) {
        thresholdsML.emplace_back(dielectroncuts.cutsMl.value[i]);
      }
      fDielectronCut.SetMLThresholds(binsML, thresholdsML);

      // static constexpr int nClassesMl = 2;
      // const std::vector<int> cutDirMl = {o2::cuts_ml::CutNot, o2::cuts_ml::CutSmaller};
      // const std::vector<std::string> labelsClasses = {"Background", "Signal"};
      // const uint32_t nBinsMl = dielectroncuts.binsMl.value.size() - 1;
      // const std::vector<std::string> labelsBins(nBinsMl, "bin");
      // double cutsMlArr[nBinsMl][nClassesMl];
      // for (uint32_t i = 0; i < nBinsMl; i++) {
      //   cutsMlArr[i][0] = 0.;
      //   cutsMlArr[i][1] = dielectroncuts.cutsMl.value[i];
      // }
      // o2::framework::LabeledArray<double> cutsMl = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      // mlResponseSingleTrack.configure(dielectroncuts.binsMl.value, cutsMl, cutDirMl, nClassesMl);
      // if (dielectroncuts.loadModelsFromCCDB) {
      //   ccdbApi.init(ccdburl);
      //   mlResponseSingleTrack.setModelPathsCCDB(dielectroncuts.onnxFileNames.value, ccdbApi, dielectroncuts.onnxPathsCCDB.value, dielectroncuts.timestampCCDB.value);
      // } else {
      //   mlResponseSingleTrack.setModelPathsLocal(dielectroncuts.onnxFileNames.value);
      // }
      // mlResponseSingleTrack.cacheInputFeaturesIndices(dielectroncuts.namesInputFeatures);
      // mlResponseSingleTrack.cacheBinningIndex(dielectroncuts.nameBinningFeature);
      // mlResponseSingleTrack.init(dielectroncuts.enableOptimizations.value);

      // fDielectronCut.SetPIDMlResponse(&mlResponseSingleTrack);
    } // end of PID ML
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // for track
    fDimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, dimuoncuts.cfg_max_pt_track);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    fDimuonCut.SetTrackPhiRange(dimuoncuts.cfg_min_phi_track, dimuoncuts.cfg_max_phi_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 20);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetChi2MFT(0.f, dimuoncuts.cfg_max_chi2mft);
    fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
    fDimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(dimuoncuts.cfg_max_relDPt_wrt_matchedMCHMID, dimuoncuts.cfg_max_DEta_wrt_matchedMCHMID, dimuoncuts.cfg_max_DPhi_wrt_matchedMCHMID); // this is relevant for global muons
    fDimuonCut.SetMFTHitMap(dimuoncuts.requireMFTHitMap, dimuoncuts.requiredMFTDisks);
  }

  template <typename TTrack>
  void fillElectronInfo(TTrack const& track)
  {
    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[track.globalIndex()];
    }

    float dca3D = dca3DinSigma(track);
    float dcaXY = dcaXYinSigma(track);
    float dcaZ = dcaZinSigma(track);
    float phiPosition = track.phi() + std::asin(-0.30282 * track.sign() * (d_bz * 0.1) * dielectroncuts.cfgRefR / (2.f * track.pt()));
    o2::math_utils::bringTo02Pi(phiPosition);

    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/positive/hs"), track.pt(), track.eta(), track.phi(), dca3D, dcaXY, dcaZ, weight);
      fRegistry.fill(HIST("Track/positive/hPhiPosition"), phiPosition);
      fRegistry.fill(HIST("Track/positive/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/positive/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/positive/hDCAxyzSigma"), dcaXY, dcaZ);
      fRegistry.fill(HIST("Track/positive/hDCAxyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/positive/hDCAzRes_Pt"), track.pt(), std::sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/positive/hDCA3dRes_Pt"), track.pt(), sigmaDca3D(track) * 1e+4);      // convert cm to um
      fRegistry.fill(HIST("Track/positive/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/positive/hNclsTPC_Pt"), track.pt(), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/positive/hNcrTPC_Pt"), track.pt(), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/positive/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/positive/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/positive/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
      if (track.hasTPC()) {
        fRegistry.fill(HIST("Track/positive/hDeltaPin"), track.tpcInnerParam(), (track.p() - track.tpcInnerParam()) / track.tpcInnerParam());
      }
      fRegistry.fill(HIST("Track/positive/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/positive/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/positive/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/positive/hChi2TOF"), track.p(), track.tofChi2());

      fRegistry.fill(HIST("Track/positive/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/positive/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/positive/hProbElBDT"), track.tpcInnerParam(), track.probElBDT());
      fRegistry.fill(HIST("Track/positive/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/positive/hMeanClusterSizeITSib"), track.p(), track.meanClusterSizeITSib() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/positive/hMeanClusterSizeITSob"), track.p(), track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      // fRegistry.fill(HIST("Track/positive/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
    } else {
      fRegistry.fill(HIST("Track/negative/hs"), track.pt(), track.eta(), track.phi(), dca3D, dcaXY, dcaZ, weight);
      fRegistry.fill(HIST("Track/negative/hPhiPosition"), phiPosition);
      fRegistry.fill(HIST("Track/negative/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/negative/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/negative/hDCAxyzSigma"), dcaXY, dcaZ);
      fRegistry.fill(HIST("Track/negative/hDCAxyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/negative/hDCAzRes_Pt"), track.pt(), std::sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/negative/hDCA3dRes_Pt"), track.pt(), sigmaDca3D(track) * 1e+4);      // convert cm to um
      fRegistry.fill(HIST("Track/negative/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/negative/hNclsTPC_Pt"), track.pt(), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/negative/hNcrTPC_Pt"), track.pt(), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/negative/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/negative/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/negative/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
      if (track.hasTPC()) {
        fRegistry.fill(HIST("Track/negative/hDeltaPin"), track.tpcInnerParam(), (track.p() - track.tpcInnerParam()) / track.tpcInnerParam());
      }
      fRegistry.fill(HIST("Track/negative/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/negative/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/negative/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/negative/hChi2TOF"), track.p(), track.tofChi2());

      fRegistry.fill(HIST("Track/negative/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/negative/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/negative/hProbElBDT"), track.tpcInnerParam(), track.probElBDT());
      fRegistry.fill(HIST("Track/negative/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/negative/hMeanClusterSizeITSib"), track.p(), track.meanClusterSizeITSib() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/negative/hMeanClusterSizeITSob"), track.p(), track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      // fRegistry.fill(HIST("Track/negative/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
    }
  }

  template <typename TTrack>
  void fillMuonInfo(TTrack const& track)
  {
    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[track.globalIndex()];
    }
    float dca_xy = fwdDcaXYinSigma(track);

    float reldpt = (track.ptMatchedMCHMID() - track.pt()) / track.pt();
    float deta = track.etaMatchedMCHMID() - track.eta();
    float dphi = track.phiMatchedMCHMID() - track.phi();
    o2::math_utils::bringToPMPi(dphi);

    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/positive/hs"), track.pt(), track.eta(), track.phi(), dca_xy, weight);
      fRegistry.fill(HIST("Track/positive/hEtaPhi_MatchMCHMID"), track.phiMatchedMCHMID(), track.etaMatchedMCHMID(), weight);
      fRegistry.fill(HIST("Track/positive/hsDelta"), track.pt(), reldpt, deta, dphi, weight);
      fRegistry.fill(HIST("Track/positive/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/positive/hTrackType"), track.trackType());
      fRegistry.fill(HIST("Track/positive/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
      fRegistry.fill(HIST("Track/positive/hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXXatDCA()), track.fwdDcaY() / std::sqrt(track.cYYatDCA()));
      fRegistry.fill(HIST("Track/positive/hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXXatDCA()) * 1e+4);
      fRegistry.fill(HIST("Track/positive/hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYYatDCA()) * 1e+4);
      fRegistry.fill(HIST("Track/positive/hDCAxyRes_Pt"), track.pt(), sigmaFwdDcaXY(track) * 1e+4);
      fRegistry.fill(HIST("Track/positive/hNclsMCH"), track.nClusters());
      fRegistry.fill(HIST("Track/positive/hNclsMFT"), track.nClustersMFT());
      fRegistry.fill(HIST("Track/positive/hPDCA"), track.rAtAbsorberEnd(), track.pDca());
      fRegistry.fill(HIST("Track/positive/hChi2"), track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) ? track.chi2() / (2.f * (track.nClusters() + track.nClustersMFT()) - 5.f) : track.chi2());
      fRegistry.fill(HIST("Track/positive/hChi2MFT"), track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) ? track.chi2MFT() / (2.f * track.nClustersMFT() - 5.f) : 0);
      fRegistry.fill(HIST("Track/positive/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
      fRegistry.fill(HIST("Track/positive/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
      fRegistry.fill(HIST("Track/positive/hMFTClusterMap"), track.mftClusterMap());
    } else {
      fRegistry.fill(HIST("Track/negative/hs"), track.pt(), track.eta(), track.phi(), dca_xy, weight);
      fRegistry.fill(HIST("Track/negative/hEtaPhi_MatchMCHMID"), track.phiMatchedMCHMID(), track.etaMatchedMCHMID(), weight);
      fRegistry.fill(HIST("Track/negative/hsDelta"), track.pt(), reldpt, deta, dphi, weight);
      fRegistry.fill(HIST("Track/negative/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/negative/hTrackType"), track.trackType());
      fRegistry.fill(HIST("Track/negative/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
      fRegistry.fill(HIST("Track/negative/hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXXatDCA()), track.fwdDcaY() / std::sqrt(track.cYYatDCA()));
      fRegistry.fill(HIST("Track/negative/hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXXatDCA()) * 1e+4);
      fRegistry.fill(HIST("Track/negative/hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYYatDCA()) * 1e+4);
      fRegistry.fill(HIST("Track/negative/hDCAxyRes_Pt"), track.pt(), sigmaFwdDcaXY(track) * 1e+4);
      fRegistry.fill(HIST("Track/negative/hNclsMCH"), track.nClusters());
      fRegistry.fill(HIST("Track/negative/hNclsMFT"), track.nClustersMFT());
      fRegistry.fill(HIST("Track/negative/hPDCA"), track.rAtAbsorberEnd(), track.pDca());
      fRegistry.fill(HIST("Track/negative/hChi2"), track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) ? track.chi2() / (2.f * (track.nClusters() + track.nClustersMFT()) - 5.f) : track.chi2());
      fRegistry.fill(HIST("Track/negative/hChi2MFT"), track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) ? track.chi2MFT() / (2.f * track.nClustersMFT() - 5.f) : 0);
      fRegistry.fill(HIST("Track/negative/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
      fRegistry.fill(HIST("Track/negative/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
      fRegistry.fill(HIST("Track/negative/hMFTClusterMap"), track.mftClusterMap());
    }
  }

  template <bool isTriggerAnalysis, typename TCollisions, typename TTracks, typename TPreslice, typename TCut>
  void runQC(TCollisions const& collisions, TTracks const& tracks, TPreslice const& perCollision, TCut const& cut)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      if constexpr (isTriggerAnalysis) {
        if (!collision.swtalias_bit(o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value))) {
          continue;
        }
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        for (const auto& track : tracks_per_coll) {
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut.template IsSelectedTrack<false>(track)) {
              continue;
            }
          } else { // cut-based
            if (!cut.template IsSelectedTrack<false>(track)) {
              continue;
            }
          }
          fillElectronInfo(track);
        } // end of track loop
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        for (const auto& track : tracks_per_coll) {
          if (!cut.template IsSelectedTrack<false>(track)) {
            continue;
          }
          if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch(track, cut, tracks)) {
            continue;
          }

          fillMuonInfo(track);
        } // end of track loop
      }
    } // end of collision loop
  }

  std::unordered_map<int, float> map_weight;
  template <bool isTriggerAnalysis, typename TCollisions, typename TTracks, typename TPreslice, typename TCut>
  void fillTrackWeightMap(TCollisions const& collisions, TTracks const& tracks, TPreslice const& perCollision, TCut const& cut)
  {
    std::vector<int> passed_trackIds;
    passed_trackIds.reserve(tracks.size());
    for (const auto& collision : collisions) {
      initCCDB(collision);
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      if constexpr (isTriggerAnalysis) {
        if (!collision.swtalias_bit(o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value))) {
          continue;
        }
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        for (const auto& track : tracks_per_coll) {
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut.template IsSelectedTrack<false>(track)) {
              continue;
            }
          } else { // cut-based
            if (!cut.template IsSelectedTrack<false>(track)) {
              continue;
            }
          }
          passed_trackIds.emplace_back(track.globalIndex());
        } // end of track loop
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        for (const auto& track : tracks_per_coll) {
          if (!cut.template IsSelectedTrack<false>(track)) {
            continue;
          }
          if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch(track, cut, tracks)) {
            continue;
          }
          passed_trackIds.emplace_back(track.globalIndex());
        } // end of track loop
      }
    } // end of collision loop

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      for (const auto& trackId : passed_trackIds) {
        auto track = tracks.rawIteratorAt(trackId);
        auto ambIds = track.ambiguousElectronsIds();
        float n = 1.f; // include myself.
        for (const auto& ambId : ambIds) {
          if (std::find(passed_trackIds.begin(), passed_trackIds.end(), ambId) != passed_trackIds.end()) {
            n += 1.f;
          }
        }
        map_weight[trackId] = 1.f / n;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      for (const auto& trackId : passed_trackIds) {
        auto track = tracks.rawIteratorAt(trackId);
        auto ambIds = track.ambiguousMuonsIds();
        float n = 1.f; // include myself.
        for (const auto& ambId : ambIds) {
          if (std::find(passed_trackIds.begin(), passed_trackIds.end(), ambId) != passed_trackIds.end()) {
            n += 1.f;
          }
        }
        map_weight[trackId] = 1.f / n;
      }
    }

    passed_trackIds.clear();
    passed_trackIds.shrink_to_fit();
  }

  SliceCache cache;
  Preslice<MyElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Filter trackFilter_electron = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter_electron = dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl;
  Filter ttcaFilter_electron = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);

  Preslice<MyMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && o2::aod::fwdtrack::pt < dimuoncuts.cfg_max_pt_track && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track && dimuoncuts.cfg_min_phi_track < o2::aod::fwdtrack::phi && o2::aod::fwdtrack::phi < dimuoncuts.cfg_max_phi_track;
  Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_numContrib = cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < cfgNumContribMax;
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processQC(FilteredMyCollisions const& collisions, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<false>(collisions, electrons, perCollision_electron, fDielectronCut);
      }
      runQC<false>(collisions, electrons, perCollision_electron, fDielectronCut);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<false>(collisions, muons, perCollision_muon, fDimuonCut);
      }
      runQC<false>(collisions, muons, perCollision_muon, fDimuonCut);
    }

    map_weight.clear();
  }
  PROCESS_SWITCH(SingleTrackQC, processQC, "run single track QC", true);

  using FilteredMyCollisionsWithSWT = soa::Filtered<MyCollisionsWithSWT>;
  void processQC_TriggeredData(FilteredMyCollisionsWithSWT const& collisions, aod::EMSWTriggerInfos const& cefpinfos, aod::EMSWTriggerATCounters const& countersAT, aod::EMSWTriggerTOICounters const& countersTOI, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<true>(collisions, electrons, perCollision_electron, fDielectronCut);
      }
      runQC<true>(collisions, electrons, perCollision_electron, fDielectronCut);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<true>(collisions, muons, perCollision_muon, fDimuonCut);
      }
      runQC<true>(collisions, muons, perCollision_muon, fDimuonCut);
    }
    map_weight.clear();

    // for nomalization
    int emswtId = o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value);
    for (const auto& counter : countersAT) {
      if (counter.isAnalyzed_bit(emswtId)) {
        fRegistry.fill(HIST("NormTrigger/hTriggerCounter"), mRunNumber, 0);
      }
    }
    for (const auto& counter : countersTOI) {
      if (counter.isAnalyzedToI_bit(emswtId)) {
        fRegistry.fill(HIST("NormTrigger/hTriggerCounter"), mRunNumber, 1);
      }
    }

    for (const auto& info : cefpinfos) {
      fRegistry.fill(HIST("NormTrigger/hInspectedTVX"), info.runNumber(), info.nInspectedTVX());
      fRegistry.fill(HIST("NormTrigger/hScalers"), info.runNumber(), info.nScalers()[emswtId]);
      fRegistry.fill(HIST("NormTrigger/hSelections"), info.runNumber(), info.nSelections()[emswtId]);
    }
  }
  PROCESS_SWITCH(SingleTrackQC, processQC_TriggeredData, "run single track QC on triggered data", false);

  void processNorm(aod::EMEventNormInfos const& collisions)
  {
    for (const auto& collision : collisions) {
      fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 1.0);
      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 2.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 3.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 4.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 5.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 6.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 7.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 8.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 9.0);
      }
      if (collision.sel8()) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 10.0);
      }
      if (std::fabs(collision.posZ()) < 10.0) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 11.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 12.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 13.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 14.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 15.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 16.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 17.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 18.0);
      }
      if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        fRegistry.fill(HIST("Event/norm/hCollisionCounter"), 19.0);
      }
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      fRegistry.fill(HIST("Event/norm/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
    } // end of collision loop
  }
  PROCESS_SWITCH(SingleTrackQC, processNorm, "process normalization info", false);

  void processBC(aod::EMBCs const& bcs)
  {
    for (const auto& bc : bcs) {
      if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fRegistry.fill(HIST("BC/hTVXCounter"), 0.f);

        if (bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
          fRegistry.fill(HIST("BC/hTVXCounter"), 1.f);
        }
        if (bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
          fRegistry.fill(HIST("BC/hTVXCounter"), 2.f);
        }
        if (rctChecker(bc)) {
          fRegistry.fill(HIST("BC/hTVXCounter"), 3.f);
        }
        if (bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
          fRegistry.fill(HIST("BC/hTVXCounter"), 4.f);
        }
        if (bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) && rctChecker(bc)) {
          fRegistry.fill(HIST("BC/hTVXCounter"), 5.f);
        }
      }
    }
  }
  PROCESS_SWITCH(SingleTrackQC, processBC, "process BC counter", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(SingleTrackQC, processDummy, "Dummy function", false);
};
#endif // PWGEM_DILEPTON_CORE_SINGLETRACKQC_H_
