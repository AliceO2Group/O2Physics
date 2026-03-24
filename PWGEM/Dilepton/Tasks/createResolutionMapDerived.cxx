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
// Analysis task to produce resolution map for electrons/muons over derived data.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <array>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;

struct createResolutionMapDerived {

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<bool> cfg_require_true_mc_collision_association{"cfg_require_true_mc_collision_association", false, "flag to require true mc collision association"};

  ConfigurableAxis ConfPtGenBins{"ConfPtGenBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00}, "gen. pT bins for output histograms"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0, 10, 30, 50, 110}, "centrality (%) bins for output histograms"};

  ConfigurableAxis ConfEtaCBGenBins{"ConfEtaCBGenBins", {30, -1.5, +1.5}, "gen. eta bins at midrapidity for output histograms"};
  ConfigurableAxis ConfEtaFWDGenBins{"ConfEtaFWDGenBins", {40, -5.5, -1.5}, "gen. eta bins at forward rapidity for output histograms"};
  ConfigurableAxis ConfPhiGenBins{"ConfPhiGenBins", {36, 0, 2.f * M_PI}, "gen. phi bins at forward rapidity for output histograms"};
  Configurable<float> cfgRefR{"cfgRefR", 0.50, "ref. radius (m) for calculating phi position"}; // 0.50 +/- 0.06 can be syst. unc.

  ConfigurableAxis ConfRelDeltaPtCBBins{"ConfRelDeltaPtCBBins", {200, -1.f, +1.f}, "rel. dpt for output histograms at midrapidity"};
  ConfigurableAxis ConfRelDeltaPtFWDBins{"ConfRelDeltaPtFWDBins", {200, -1.f, +1.f}, "rel. dpt for output histograms at fwd rapidity"};

  ConfigurableAxis ConfDeltaEtaCBBins{"ConfDeltaEtaCBBins", {200, -0.5f, +0.5f}, "deta bins for output histograms at midrapidity"};
  ConfigurableAxis ConfDeltaEtaFWDBins{"ConfDeltaEtaFWDBins", {200, -0.5f, +0.5f}, "deta bins for output histograms at fwd rapidity"};
  ConfigurableAxis ConfDeltaPhiBins{"ConfDeltaPhiBins", {200, -0.5f, +0.5f}, "dphi bins for output histograms"};

  Configurable<bool> cfgFillTHnSparse{"cfgFillTHnSparse", true, "fill THnSparse for output"};
  Configurable<bool> cfgFillTH2{"cfgFillTH2", false, "fill TH2 for output"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
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
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT_hadronPID, CBT_muon_glo, CBT_muon]. see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb/OO"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

    Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
    Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
    Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
  } eventcuts;

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "electroncut_group";

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.05, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -2, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +2, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<bool> cfg_mirror_phi_track{"cfg_mirror_phi_track", false, "mirror the phi cut around Pi, min and max Phi should be in 0-Pi"};
    Configurable<bool> cfg_reject_phi_track{"cfg_reject_phi_track", false, "reject the phi interval"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2 TOF"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_min_its_cluster_size{"cfg_min_its_cluster_size", 0.f, "min ITS cluster size"};
    Configurable<float> cfg_max_its_cluster_size{"cfg_max_its_cluster_size", 16.f, "max ITS cluster size"};
    Configurable<float> cfgRefR{"cfgRefR", 0.50, "ref. radius (m) for calculating phi position"}; // 0.50 +/- 0.06 can be syst. unc.
    Configurable<float> cfg_min_phiposition_track{"cfg_min_phiposition_track", 0.f, "min phi position for single track at certain radius"};
    Configurable<float> cfg_max_phiposition_track{"cfg_max_phiposition_track", 6.3, "max phi position for single track at certain radius"};
    Configurable<bool> acceptOnlyCorrectMatch{"acceptOnlyCorrectMatch", false, "flag to accept only correct match between ITS and TPC"}; // this is only for MC study, as we don't know correct match in data.
    Configurable<bool> acceptOnlyWrongMatch{"acceptOnlyWrongMatch", false, "flag to accept only wrong match between ITS and TPC"};       // this is only for MC study, as we don't know correct match in data.

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

    Configurable<bool> checkPIDForTracking{"checkPIDForTracking", false, "check PID for tracking"};
    Configurable<uint32_t> PartIdentifier{"PartIdentifier", 2, "Particle identifier for selected particle; 0: electron, 1: muon, 2: pion, 3: kaon, 4: proton, 5: deuteron, 6: triton, 7: helium3, 8: alpha"};

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
  } electroncuts;

  DimuonCut fDimuonCut;
  struct : ConfigurableGroup {
    std::string prefix = "muoncut_group";

    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 0, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -10, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "max phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+6, "max chi2/ndf"};
    Configurable<float> cfg_max_chi2mft{"cfg_max_chi2mft", 1e+6, "max chi2/ndf"};
    // Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 40, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_border_pt_for_chi2mchmft{"cfg_border_pt_for_chi2mchmft", 0, "border pt for different max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mftmch_lowPt{"cfg_max_matching_chi2_mftmch_lowPt", 8, "max chi2 for MFT-MCH matching for low pT"};
    Configurable<float> cfg_max_matching_chi2_mftmch_highPt{"cfg_max_matching_chi2_mftmch_highPt", 40, "max chi2 for MFT-MCH matching for high pT"};
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
    Configurable<bool> acceptOnlyCorrectMatch{"acceptOnlyCorrectMatch", false, "flag to accept only correct match between MFT and MCH-MID"}; // this is only for MC study, as we don't know correct match in data.
    Configurable<bool> acceptOnlyWrongMatch{"acceptOnlyWrongMatch", false, "flag to accept only wrong match between MFT and MCH-MID"};       // this is only for MC study, as we don't know correct match in data.
  } muoncuts;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  ~createResolutionMapDerived() {}

  void init(o2::framework::InitContext&)
  {
    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);

    DefineEMEventCut();

    if (doprocessElectron) {
      DefineDielectronCut();
    }
    if (doprocessMuon) {
      DefineDimuonCut();
    }
    addHistograms();
  }

  void addHistograms()
  {
    // registry.add("Event/Electron/hImpPar_Centrality", "true imapact parameter vs. estimated centrality;impact parameter (fm);centrality (%)", kTH2F, {{200, 0, 20}, {110, 0, 110}}, true);
    // registry.add("Event/Electron/hImpPar_Centrality", "true imapact parameter vs. estimated centrality;impact parameter (fm);centrality (%)", kTH2F, {{200, 0, 20}, {110, 0, 110}}, true);

    if (doprocessGen) {
      registry.add("Event/hGenID", "generator ID;generator ID;Number of mc collisions", kTH1F, {{7, -1.5, 5.5}}, true);
    }

    const AxisSpec axis_cent{ConfCentBins, "centrality (%)"};
    const AxisSpec axis_pt_gen{ConfPtGenBins, "p_{T,l}^{gen} (GeV/c)"};
    const AxisSpec axis_eta_cb_gen{ConfEtaCBGenBins, "#eta_{l}^{gen}"};
    const AxisSpec axis_eta_fwd_gen{ConfEtaFWDGenBins, "#eta_{l}^{gen}"};
    const AxisSpec axis_phi_gen{ConfPhiGenBins, "#varphi_{l}^{gen} (rad.)"};
    const AxisSpec axis_dpt_cb{ConfRelDeltaPtCBBins, "(p_{T,l}^{gen} - p_{T,l}^{rec})/p_{T,l}^{gen}"};
    const AxisSpec axis_dpt_fwd{ConfRelDeltaPtFWDBins, "(p_{T,l}^{gen} - p_{T,l}^{rec})/p_{T,l}^{gen}"};
    const AxisSpec axis_deta_cb{ConfDeltaEtaCBBins, "#eta_{l}^{gen} - #eta_{l}^{rec}"};
    const AxisSpec axis_deta_fwd{ConfDeltaEtaFWDBins, "#eta_{l}^{gen} - #eta_{l}^{rec}"};
    const AxisSpec axis_dphi{ConfDeltaPhiBins, "#varphi_{l}^{gen} - #varphi_{l}^{rec} (rad.)"};
    const AxisSpec axis_charge_gen{3, -1.5, +1.5, "true sign"};

    if (doprocessElectron) {
      registry.add("Electron/hPIDForTracking", "PID for trackng", kTH1F, {{9, -0.5, 8.5}}, false); // see numbering in O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h
      if (cfgFillTH2) {
        registry.add("Electron/hPt", "rec. p_{T,e};p_{T,e} (GeV/c)", kTH1F, {{1000, 0, 10}}, false);
        registry.add("Electron/hEtaPhi", "rec. #eta vs. #varphi;#varphi_{e} (rad.);#eta_{e}", kTH2F, {{90, 0, 2 * M_PI}, {40, -2, +2}}, false);
        registry.add("Electron/Ptgen_RelDeltaPt", "resolution", kTH2F, {{axis_pt_gen}, {axis_dpt_cb}}, true);
        registry.add("Electron/Ptgen_DeltaEta", "resolution", kTH2F, {{axis_pt_gen}, {axis_deta_cb}}, true);
        registry.add("Electron/Ptgen_DeltaPhi_Pos", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
        registry.add("Electron/Ptgen_DeltaPhi_Neg", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
      }
      if (cfgFillTHnSparse) {
        registry.add("Electron/hs_reso", "8D resolution", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_cb_gen, axis_phi_gen, axis_charge_gen, axis_dpt_cb, axis_deta_cb, axis_dphi}, true);
      }
    }

    if (doprocessMuon) {
      if (muoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
        if (cfgFillTH2) {
          registry.add("GlobalMuon/hPt", "rec. p_{T,#mu};p_{T,#mu} (GeV/c)", kTH1F, {{1000, 0, 10}}, false);
          registry.add("GlobalMuon/hEtaPhi", "rec. #eta vs. #varphi;#varphi_{#mu} (rad.);#eta_{#mu}", kTH2F, {{90, 0, 2 * M_PI}, {60, -6, 0}}, false);
          registry.add("GlobalMuon/Ptgen_RelDeltaPt", "resolution", kTH2F, {{axis_pt_gen}, {axis_dpt_fwd}}, true);
          registry.add("GlobalMuon/Ptgen_DeltaEta", "resolution", kTH2F, {{axis_pt_gen}, {axis_deta_fwd}}, true);
          registry.add("GlobalMuon/Ptgen_DeltaPhi_Pos", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
          registry.add("GlobalMuon/Ptgen_DeltaPhi_Neg", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
        }
        if (cfgFillTHnSparse) {
          registry.add("GlobalMuon/hs_reso", "8D resolution", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt_fwd, axis_deta_fwd, axis_dphi}, true);
        }
      } else if (muoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
        if (cfgFillTH2) {
          registry.add("StandaloneMuon/hPt", "rec. p_{T,#mu};p_{T,#mu} (GeV/c)", kTH1F, {{1000, 0, 10}}, false);
          registry.add("StandaloneMuon/hEtaPhi", "rec. #eta vs. #varphi;#varphi_{#mu} (rad.);#eta_{#mu}", kTH2F, {{90, 0, 2 * M_PI}, {60, -6, 0}}, false);
          registry.add("StandaloneMuon/Ptgen_RelDeltaPt", "resolution", kTH2F, {{axis_pt_gen}, {axis_dpt_fwd}}, true);
          registry.add("StandaloneMuon/Ptgen_DeltaEta", "resolution", kTH2F, {{axis_pt_gen}, {axis_deta_fwd}}, true);
          registry.add("StandaloneMuon/Ptgen_DeltaPhi_Pos", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
          registry.add("StandaloneMuon/Ptgen_DeltaPhi_Neg", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
        }
        if (cfgFillTHnSparse) {
          registry.add("StandaloneMuon/hs_reso", "8D resolution", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt_fwd, axis_deta_fwd, axis_dphi}, true);
        }
      }
    }
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

  void DefineDielectronCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for track
    fDielectronCut.SetTrackPtRange(electroncuts.cfg_min_pt_track, electroncuts.cfg_max_pt_track);
    fDielectronCut.SetTrackEtaRange(electroncuts.cfg_min_eta_track, +electroncuts.cfg_max_eta_track);
    fDielectronCut.SetTrackPhiRange(electroncuts.cfg_min_phi_track, electroncuts.cfg_max_phi_track, electroncuts.cfg_mirror_phi_track, electroncuts.cfg_reject_phi_track);
    fDielectronCut.SetMinNClustersTPC(electroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(electroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetMaxFracSharedClustersTPC(electroncuts.cfg_max_frac_shared_clusters_tpc);
    fDielectronCut.SetChi2PerClusterTPC(0.0, electroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, electroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(electroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITS(electroncuts.cfg_min_its_cluster_size, electroncuts.cfg_max_its_cluster_size);
    fDielectronCut.SetTrackMaxDcaXY(electroncuts.cfg_max_dcaxy);
    fDielectronCut.SetTrackMaxDcaZ(electroncuts.cfg_max_dcaz);
    fDielectronCut.RequireITSibAny(electroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(electroncuts.cfg_require_itsib_1st);
    fDielectronCut.SetChi2TOF(0.0, electroncuts.cfg_max_chi2tof);
    // fDielectronCut.SetRelDiffPin(electroncuts.cfg_min_rel_diff_pin, electroncuts.cfg_max_rel_diff_pin);
    fDielectronCut.EnableTTCA(electroncuts.enableTTCA);

    // for eID
    fDielectronCut.SetPIDScheme(electroncuts.cfg_pid_scheme);
    fDielectronCut.SetTPCNsigmaElRange(electroncuts.cfg_min_TPCNsigmaEl, electroncuts.cfg_max_TPCNsigmaEl);
    // fDielectronCut.SetTPCNsigmaMuRange(electroncuts.cfg_min_TPCNsigmaMu, electroncuts.cfg_max_TPCNsigmaMu);
    fDielectronCut.SetTPCNsigmaPiRange(electroncuts.cfg_min_TPCNsigmaPi, electroncuts.cfg_max_TPCNsigmaPi);
    fDielectronCut.SetTPCNsigmaKaRange(electroncuts.cfg_min_TPCNsigmaKa, electroncuts.cfg_max_TPCNsigmaKa);
    fDielectronCut.SetTPCNsigmaPrRange(electroncuts.cfg_min_TPCNsigmaPr, electroncuts.cfg_max_TPCNsigmaPr);
    fDielectronCut.SetTOFNsigmaElRange(electroncuts.cfg_min_TOFNsigmaEl, electroncuts.cfg_max_TOFNsigmaEl);
    fDielectronCut.SetPinRangeForPionRejectionTPC(electroncuts.cfg_min_pin_pirejTPC, electroncuts.cfg_max_pin_pirejTPC);

    if (electroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      std::vector<float> binsML{};
      binsML.reserve(electroncuts.binsMl.value.size());
      for (size_t i = 0; i < electroncuts.binsMl.value.size(); i++) {
        binsML.emplace_back(electroncuts.binsMl.value[i]);
      }
      std::vector<float> thresholdsML{};
      thresholdsML.reserve(electroncuts.cutsMl.value.size());
      for (size_t i = 0; i < electroncuts.cutsMl.value.size(); i++) {
        thresholdsML.emplace_back(electroncuts.cutsMl.value[i]);
      }
      fDielectronCut.SetMLThresholds(binsML, thresholdsML);

      // static constexpr int nClassesMl = 2;
      // const std::vector<int> cutDirMl = {o2::cuts_ml::CutNot, o2::cuts_ml::CutSmaller};
      // const std::vector<std::string> labelsClasses = {"Background", "Signal"};
      // const uint32_t nBinsMl = electroncuts.binsMl.value.size() - 1;
      // const std::vector<std::string> labelsBins(nBinsMl, "bin");
      // double cutsMlArr[nBinsMl][nClassesMl];
      // for (uint32_t i = 0; i < nBinsMl; i++) {
      //   cutsMlArr[i][0] = 0.;
      //   cutsMlArr[i][1] = electroncuts.cutsMl.value[i];
      // }
      // o2::framework::LabeledArray<double> cutsMl = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      // mlResponseSingleTrack.configure(electroncuts.binsMl.value, cutsMl, cutDirMl, nClassesMl);
      // if (electroncuts.loadModelsFromCCDB) {
      //   ccdbApi.init(ccdburl);
      //   mlResponseSingleTrack.setModelPathsCCDB(electroncuts.onnxFileNames.value, ccdbApi, electroncuts.onnxPathsCCDB.value, electroncuts.timestampCCDB.value);
      // } else {
      //   mlResponseSingleTrack.setModelPathsLocal(electroncuts.onnxFileNames.value);
      // }
      // mlResponseSingleTrack.cacheInputFeaturesIndices(electroncuts.namesInputFeatures);
      // mlResponseSingleTrack.cacheBinningIndex(electroncuts.nameBinningFeature);
      // mlResponseSingleTrack.init(electroncuts.enableOptimizations.value);

    } // end of PID ML
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // for track
    fDimuonCut.SetTrackType(muoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(muoncuts.cfg_min_pt_track, muoncuts.cfg_max_pt_track);
    fDimuonCut.SetTrackEtaRange(muoncuts.cfg_min_eta_track, muoncuts.cfg_max_eta_track);
    fDimuonCut.SetTrackPhiRange(muoncuts.cfg_min_phi_track, muoncuts.cfg_max_phi_track);
    fDimuonCut.SetNClustersMFT(muoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(muoncuts.cfg_min_ncluster_mch, 20);
    fDimuonCut.SetChi2(0.f, muoncuts.cfg_max_chi2);
    fDimuonCut.SetChi2MFT(0.f, muoncuts.cfg_max_chi2mft);
    // fDimuonCut.SetMatchingChi2MCHMFT(0.f, muoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMaxMatchingChi2MCHMFTPtDep([&](float pt) { return (pt < muoncuts.cfg_border_pt_for_chi2mchmft ? muoncuts.cfg_max_matching_chi2_mftmch_lowPt : muoncuts.cfg_max_matching_chi2_mftmch_highPt); });
    fDimuonCut.SetMatchingChi2MCHMID(0.f, muoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, muoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(muoncuts.cfg_min_rabs, muoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
    fDimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(muoncuts.cfg_max_relDPt_wrt_matchedMCHMID, muoncuts.cfg_max_DEta_wrt_matchedMCHMID, muoncuts.cfg_max_DPhi_wrt_matchedMCHMID); // this is relevant for global muons
    fDimuonCut.SetMFTHitMap(muoncuts.requireMFTHitMap, muoncuts.requiredMFTDisks);
    fDimuonCut.EnableTTCA(muoncuts.enableTTCA);
  }

  template <typename TCollisions, typename TTracks, typename TPerCollision, typename TCut>
  void create(TCollisions const& collisions, TTracks const& tracks, TPerCollision const& perCollision, TCut const& cut)
  {
    for (const auto& collision : collisions) {
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      float centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()}[eventcuts.cfgCentEstimator];
      if (centrality < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centrality) {
        continue;
      }

      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        return;
      }

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& track : tracks_per_coll) {
        if (!cut.IsSelectedTrack(track)) {
          continue;
        }

        auto mcParticle = track.template emmcparticle_as<aod::EMMCParticles>();
        auto mcCollision = mcParticle.template emmcevent_as<aod::EMMCEvents>();
        if (cfgEventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }

        if (cfg_require_true_mc_collision_association && mcParticle.emmceventId() != collision.emmceventId()) {
          return;
        }

        if constexpr (std::is_same_v<std::decay_t<TTracks>, MyMCElectrons>) {
          if (std::abs(mcParticle.pdgCode()) != 11) {
            continue;
          }

          if (electroncuts.acceptOnlyCorrectMatch && o2::aod::pwgem::dilepton::utils::mcutil::hasFakeMatchITSTPC(track)) {
            continue;
          }
          if (electroncuts.acceptOnlyWrongMatch && !o2::aod::pwgem::dilepton::utils::mcutil::hasFakeMatchITSTPC(track)) {
            continue;
          }

          registry.fill(HIST("Electron/hPIDForTracking"), track.pidForTracking());
          if (electroncuts.checkPIDForTracking && track.pidForTracking() != static_cast<uint32_t>(electroncuts.PartIdentifier)) {
            continue;
          }

          if (cfgFillTHnSparse) {
            registry.fill(HIST("Electron/hs_reso"), centrality, mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), -mcParticle.pdgCode() / 11, (mcParticle.pt() - track.pt()) / mcParticle.pt(), mcParticle.eta() - track.eta(), mcParticle.phi() - track.phi());
          }
          if (cfgFillTH2) {
            registry.fill(HIST("Electron/hPt"), track.pt());
            registry.fill(HIST("Electron/hEtaPhi"), track.phi(), track.eta());
            registry.fill(HIST("Electron/Ptgen_RelDeltaPt"), mcParticle.pt(), (mcParticle.pt() - track.pt()) / mcParticle.pt());
            registry.fill(HIST("Electron/Ptgen_DeltaEta"), mcParticle.pt(), mcParticle.eta() - track.eta());
            if (mcParticle.pdgCode() == -11) { // positron
              registry.fill(HIST("Electron/Ptgen_DeltaPhi_Pos"), mcParticle.pt(), mcParticle.phi() - track.phi());
            } else if (mcParticle.pdgCode() == 11) { // electron
              registry.fill(HIST("Electron/Ptgen_DeltaPhi_Neg"), mcParticle.pt(), mcParticle.phi() - track.phi());
            }
          }
        } else if constexpr (std::is_same_v<std::decay_t<TTracks>, MyMCMuons>) {
          if (std::abs(mcParticle.pdgCode()) != 13) {
            continue;
          }

          if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
            if (muoncuts.acceptOnlyCorrectMatch) {
              if (track.emmcparticleId() != track.emmftmcparticleId()) {
                continue;
              }
            }
            if (muoncuts.acceptOnlyWrongMatch) { // reject correctly matched MFT-MCH-MID for bkg estimation
              if (track.emmcparticleId() == track.emmftmcparticleId()) {
                continue;
              }
            }

            if (cfgFillTHnSparse) {
              registry.fill(HIST("GlobalMuon/hs_reso"), centrality, mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), -mcParticle.pdgCode() / 13, (mcParticle.pt() - track.pt()) / mcParticle.pt(), mcParticle.eta() - track.eta(), mcParticle.phi() - track.phi());
            }

            if (cfgFillTH2) {
              registry.fill(HIST("GlobalMuon/hPt"), track.pt());
              registry.fill(HIST("GlobalMuon/hEtaPhi"), track.phi(), track.eta());
              registry.fill(HIST("GlobalMuon/Ptgen_RelDeltaPt"), mcParticle.pt(), (mcParticle.pt() - track.pt()) / mcParticle.pt());
              registry.fill(HIST("GlobalMuon/Ptgen_DeltaEta"), mcParticle.pt(), mcParticle.eta() - track.eta());
              if (mcParticle.pdgCode() == -13) { // positive muon
                registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Pos"), mcParticle.pt(), mcParticle.phi() - track.phi());
              } else if (mcParticle.pdgCode() == 13) { // negative muon
                registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Neg"), mcParticle.pt(), mcParticle.phi() - track.phi());
              }
            }
          } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
            if (cfgFillTHnSparse) {
              registry.fill(HIST("StandaloneMuon/hs_reso"), centrality, mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), -mcParticle.pdgCode() / 13, (mcParticle.pt() - track.pt()) / mcParticle.pt(), mcParticle.eta() - track.eta(), mcParticle.phi() - track.phi());
            }
            if (cfgFillTH2) {
              registry.fill(HIST("StandaloneMuon/hPt"), track.pt());
              registry.fill(HIST("StandaloneMuon/hEtaPhi"), track.phi(), track.eta());
              registry.fill(HIST("StandaloneMuon/Ptgen_RelDeltaPt"), mcParticle.pt(), (mcParticle.pt() - track.pt()) / mcParticle.pt());
              registry.fill(HIST("StandaloneMuon/Ptgen_DeltaEta"), mcParticle.pt(), mcParticle.eta() - track.eta());
              if (mcParticle.pdgCode() == -13) { // positive muon
                registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Pos"), mcParticle.pt(), mcParticle.phi() - track.phi());
              } else if (mcParticle.pdgCode() == 13) { // negative muon
                registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Neg"), mcParticle.pt(), mcParticle.phi() - track.phi());
              }
            }
          }
        }
      } // end of track loop per collision
    } // end of collisions
  }

  using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
  using MyCollision = MyCollisions::iterator;

  using MyMCElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit, aod::EMPrimaryElectronMCLabels>;
  using MyMCElectron = MyMCElectrons::iterator;

  using MyMCMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMGlobalMuonSelfIds, aod::EMPrimaryMuonMCLabels, aod::EMMFTMCLabels>;
  using MyMCMuon = MyMCMuons::iterator;

  SliceCache cache;
  Preslice<MyMCElectrons> perCollision_mid = aod::emprimaryelectron::emeventId;
  Preslice<MyMCMuons> perCollision_fwd = aod::emprimarymuon::emeventId;

  Filter collisionFilter_centrality = (eventcuts.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventcuts.cfgCentMax);
  Filter collisionFilter_numContrib = eventcuts.cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < eventcuts.cfgNumContribMax;
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using filteredMyCollisions = soa::Filtered<MyCollisions>;

  template <typename TTracks>
  void processReso(filteredMyCollisions const& collisions, TTracks const& tracks, aod::EMMCEvents const&, aod::EMMCParticles const&)
  {
    if constexpr (std::is_same_v<std::decay_t<TTracks>, MyMCElectrons>) {
      create(collisions, tracks, perCollision_mid, fDielectronCut);
    } else if constexpr (std::is_same_v<std::decay_t<TTracks>, MyMCMuons>) {
      create(collisions, tracks, perCollision_fwd, fDimuonCut);
    }
  }

  PROCESS_SWITCH_FULL(createResolutionMapDerived, processReso<MyMCElectrons>, processElectron, "create resolution map for electrons at mid rapidity", false);
  PROCESS_SWITCH_FULL(createResolutionMapDerived, processReso<MyMCMuons>, processMuon, "create resolution map for global muons at fwd rapidity", false); // gl or sa can be selected in subwagons.

  void processGen(aod::EMMCEvents const& mcCollisions)
  {
    for (const auto& mccollision : mcCollisions) {
      registry.fill(HIST("Event/hGenID"), mccollision.getSubGeneratorId());
    }
  }
  PROCESS_SWITCH(createResolutionMapDerived, processGen, "process generated info", true);

  void processDummy(aod::EMMCEvents const&) {}
  PROCESS_SWITCH(createResolutionMapDerived, processDummy, "process dummy", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<createResolutionMapDerived>(cfgc, TaskName{"create-resolution-map-derived"})};
}
