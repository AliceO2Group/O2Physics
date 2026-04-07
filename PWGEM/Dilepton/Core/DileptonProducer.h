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
// This code runs loop over leptons.
//    Please write to: daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_CORE_DILEPTONPRODUCER_H_
#define PWGEM_DILEPTON_CORE_DILEPTONPRODUCER_H_

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <sys/types.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

using MyCollisions = o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMEventsQvec2>;
using MyCollision = MyCollisions::iterator;

using MyElectrons = o2::soa::Join<o2::aod::EMPrimaryElectrons, o2::aod::EMPrimaryElectronEMEventIds, o2::aod::EMAmbiguousElectronSelfIds, o2::aod::EMPrimaryElectronsPrefilterBit, o2::aod::EMPrimaryElectronsPrefilterBitDerived>;
using MyElectron = MyElectrons::iterator;
using FilteredMyElectrons = o2::soa::Filtered<MyElectrons>;
using FilteredMyElectron = FilteredMyElectrons::iterator;

using MyMuons = o2::soa::Join<o2::aod::EMPrimaryMuons, o2::aod::EMPrimaryMuonEMEventIds, o2::aod::EMAmbiguousMuonSelfIds, o2::aod::EMGlobalMuonSelfIds, o2::aod::EMPrimaryMuonsPrefilterBitDerived>;
using MyMuon = MyMuons::iterator;
using FilteredMyMuons = o2::soa::Filtered<MyMuons>;
using FilteredMyMuon = FilteredMyMuons::iterator;

template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename... Types>
struct DileptonProducer {
  o2::framework::Produces<o2::aod::EMThinEvents> eventTable;
  o2::framework::Produces<o2::aod::EMThinEventNormInfos> normTable;
  o2::framework::Produces<o2::aod::EMDileptons> dileptonTable;

  // Configurables
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  o2::framework::Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  o2::framework::Configurable<int> cfgEP2Estimator_for_Mix{"cfgEP2Estimator_for_Mix", 3, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  o2::framework::Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  o2::framework::Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  o2::framework::Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  o2::framework::Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};
  o2::framework::Configurable<uint> cfgDCAType{"cfgDCAType", 0, "type of DCA for output. 0:3D, 1:XY, 2:Z, else:3D"};
  o2::framework::Configurable<bool> cfgStoreULS{"cfgStoreULS", true, "flag to store ULS pairs"};
  o2::framework::Configurable<bool> cfgStoreLS{"cfgStoreLS", true, "flag to store LS pairs"};

  EMEventCut fEMEventCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "eventcut_group";
    o2::framework::Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    o2::framework::Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    o2::framework::Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    o2::framework::Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    o2::framework::Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    o2::framework::Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    o2::framework::Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    o2::framework::Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    o2::framework::Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    o2::framework::Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    o2::framework::Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    o2::framework::Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    o2::framework::Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    o2::framework::Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
    o2::framework::Configurable<bool> cfgRequireGoodITSLayer3{"cfgRequireGoodITSLayer3", false, "number of inactive chips on ITS layer 3 are below threshold "};
    o2::framework::Configurable<bool> cfgRequireGoodITSLayer0123{"cfgRequireGoodITSLayer0123", false, "number of inactive chips on ITS layers 0-3 are below threshold "};
    o2::framework::Configurable<bool> cfgRequireGoodITSLayersAll{"cfgRequireGoodITSLayersAll", false, "number of inactive chips on all ITS layers are below threshold "};
    // for RCT
    o2::framework::Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
    o2::framework::Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    o2::framework::Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
    o2::framework::Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

    o2::framework::Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    o2::framework::Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    o2::framework::Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
    o2::framework::Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
  } eventcuts;

  DielectronCut fDielectronCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "dielectroncut_group";
    o2::framework::Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    o2::framework::Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    o2::framework::Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pT"};
    o2::framework::Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pT"};
    o2::framework::Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -0.8, "min pair rapidity"};
    o2::framework::Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", +0.8, "max pair rapidity"};
    o2::framework::Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    o2::framework::Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    o2::framework::Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    o2::framework::Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    o2::framework::Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};
    o2::framework::Configurable<float> cfg_min_phiv{"cfg_min_phiv", 0.0, "min phiv (constant)"};
    o2::framework::Configurable<float> cfg_max_phiv{"cfg_max_phiv", 3.2, "max phiv (constant)"};
    o2::framework::Configurable<bool> cfg_apply_detadphi{"cfg_apply_detadphi", false, "flag to apply deta-dphi elliptic cut at PV"};
    o2::framework::Configurable<float> cfg_min_deta{"cfg_min_deta", 0.02, "min deta between 2 electrons (elliptic cut)"};
    o2::framework::Configurable<float> cfg_min_dphi{"cfg_min_dphi", 0.2, "min dphi between 2 electrons (elliptic cut)"};
    o2::framework::Configurable<float> cfg_min_opang{"cfg_min_opang", 0.0, "min opening angle"};
    o2::framework::Configurable<float> cfg_max_opang{"cfg_max_opang", 6.4, "max opening angle"};
    // o2::framework::Configurable<bool> cfg_require_diff_sides{"cfg_require_diff_sides", false, "flag to require 2 tracks are from different sides."};

    o2::framework::Configurable<bool> cfg_apply_cuts_from_prefilter{"cfg_apply_cuts_from_prefilter", false, "flag to apply prefilter set when producing derived data"};
    o2::framework::Configurable<uint16_t> cfg_prefilter_bits{"cfg_prefilter_bits", 0, "prefilter bits [kNone : 0, kElFromPC : 1, kElFromPi0_20MeV : 2, kElFromPi0_40MeV : 4, kElFromPi0_60MeV : 8, kElFromPi0_80MeV : 16, kElFromPi0_100MeV : 32, kElFromPi0_120MeV : 64, kElFromPi0_140MeV : 128] Please consider logical-OR among them."}; // see PairUtilities.h

    o2::framework::Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply pair cut same as prefilter set in derived data"};
    o2::framework::Configurable<uint16_t> cfg_prefilter_bits_derived{"cfg_prefilter_bits_derived", 0, "prefilter bits [kNone : 0, kMee : 1, kPhiV : 2, kSplitOrMergedTrackLS : 4, kSplitOrMergedTrackULS : 8] Please consider logical-OR among them."}; // see PairUtilities.h

    o2::framework::Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    o2::framework::Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    o2::framework::Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "min eta for single track"};
    o2::framework::Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
    o2::framework::Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    o2::framework::Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    o2::framework::Configurable<bool> cfg_mirror_phi_track{"cfg_mirror_phi_track", false, "mirror the phi cut around Pi, min and max Phi should be in 0-Pi"};
    o2::framework::Configurable<bool> cfg_reject_phi_track{"cfg_reject_phi_track", false, "reject the phi interval"};
    o2::framework::Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    o2::framework::Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    o2::framework::Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    o2::framework::Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2 TOF"};
    o2::framework::Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.2, "max dca XY for single track in cm"};
    o2::framework::Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.2, "max dca Z for single track in cm"};
    o2::framework::Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    o2::framework::Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    o2::framework::Configurable<float> cfg_min_its_cluster_size{"cfg_min_its_cluster_size", 0.f, "min ITS cluster size"};
    o2::framework::Configurable<float> cfg_max_its_cluster_size{"cfg_max_its_cluster_size", 16.f, "max ITS cluster size"};
    // o2::framework::Configurable<float> cfg_min_rel_diff_pin{"cfg_min_rel_diff_pin", -1e+10, "min rel. diff. between pin and ppv"};
    // o2::framework::Configurable<float> cfg_max_rel_diff_pin{"cfg_max_rel_diff_pin", +1e+10, "max rel. diff. between pin and ppv"};
    o2::framework::Configurable<float> cfgRefR{"cfgRefR", 0.50, "ref. radius (m) for calculating phi position"}; // 0.50 +/- 0.06 can be syst. unc.
    o2::framework::Configurable<float> cfg_min_phiposition_track{"cfg_min_phiposition_track", 0.f, "min phi position for single track at certain radius"};
    o2::framework::Configurable<float> cfg_max_phiposition_track{"cfg_max_phiposition_track", 6.3, "max phi position for single track at certain radius"};

    o2::framework::Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3, kTOFif : 4, kPIDML : 5, kTPChadrejORTOFreq_woTOFif : 6]"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    // o2::framework::Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    // o2::framework::Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
    o2::framework::Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    o2::framework::Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
    o2::framework::Configurable<float> cfg_min_pin_pirejTPC{"cfg_min_pin_pirejTPC", 0.f, "min. pin for pion rejection in TPC"};
    o2::framework::Configurable<float> cfg_max_pin_pirejTPC{"cfg_max_pin_pirejTPC", 1e+10, "max. pin for pion rejection in TPC"};
    o2::framework::Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};

    // configuration for PID ML
    o2::framework::Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
    o2::framework::Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
    o2::framework::Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{0.1, 0.15, 0.2, 0.25, 0.4, 0.8, 1.6, 2.0, 20.f}, "Bin limits for ML application"};
    o2::framework::Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.98, 0.98, 0.9, 0.9, 0.95, 0.95, 0.8, 0.8}, "ML cuts per bin"};
    o2::framework::Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature"}, "Names of ML model input features"};
    o2::framework::Configurable<std::string> nameBinningFeature{"nameBinningFeature", "pt", "Names of ML model binning feature"};
    o2::framework::Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    o2::framework::Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    o2::framework::Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dielectroncuts;

  DimuonCut fDimuonCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "dimuoncut_group";
    o2::framework::Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    o2::framework::Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    o2::framework::Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pt"};
    o2::framework::Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pt"};
    o2::framework::Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -4.0, "min pair rapidity"};
    o2::framework::Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", -2.5, "max pair rapidity"};
    o2::framework::Configurable<float> cfg_min_pair_dcaxy{"cfg_min_pair_dcaxy", 0.0, "min pair dca3d in sigma"};
    o2::framework::Configurable<float> cfg_max_pair_dcaxy{"cfg_max_pair_dcaxy", 1e+10, "max pair dca3d in sigma"};
    o2::framework::Configurable<bool> cfg_apply_detadphi{"cfg_apply_detadphi", false, "flag to apply deta-dphi elliptic cut"};
    o2::framework::Configurable<float> cfg_min_deta{"cfg_min_deta", 0.02, "min deta between 2 muons (elliptic cut)"};
    o2::framework::Configurable<float> cfg_min_dphi{"cfg_min_dphi", 0.02, "min dphi between 2 muons (elliptic cut)"};

    o2::framework::Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter set in derived data"};
    o2::framework::Configurable<uint16_t> cfg_prefilter_bits_derived{"cfg_prefilter_bits_derived", 0, "prefilter bits [kNone : 0, kSplitOrMergedTrackLS : 4, kSplitOrMergedTrackULS : 8] Please consider logical-OR among them."}; // see PairUtilities.h

    o2::framework::Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    o2::framework::Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    o2::framework::Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    o2::framework::Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    o2::framework::Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    o2::framework::Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    o2::framework::Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    o2::framework::Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    o2::framework::Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    o2::framework::Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+6, "max chi2/ndf"};
    o2::framework::Configurable<float> cfg_max_chi2mft{"cfg_max_chi2mft", 1e+6, "max chi2/ndf"};
    // o2::framework::Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 40, "max chi2 for MFT-MCH matching"};
    o2::framework::Configurable<float> cfg_border_pt_for_chi2mchmft{"cfg_border_pt_for_chi2mchmft", 0, "border pt for different max chi2 for MFT-MCH matching"};
    o2::framework::Configurable<float> cfg_max_matching_chi2_mftmch_lowPt{"cfg_max_matching_chi2_mftmch_lowPt", 8, "max chi2 for MFT-MCH matching for low pT"};
    o2::framework::Configurable<float> cfg_max_matching_chi2_mftmch_highPt{"cfg_max_matching_chi2_mftmch_highPt", 40, "max chi2 for MFT-MCH matching for high pT"};
    o2::framework::Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    o2::framework::Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    o2::framework::Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    o2::framework::Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    o2::framework::Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
    o2::framework::Configurable<float> cfg_max_relDPt_wrt_matchedMCHMID{"cfg_max_relDPt_wrt_matchedMCHMID", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<float> cfg_max_DEta_wrt_matchedMCHMID{"cfg_max_DEta_wrt_matchedMCHMID", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<float> cfg_max_DPhi_wrt_matchedMCHMID{"cfg_max_DPhi_wrt_matchedMCHMID", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "flag to apply MFT hit map"};
    o2::framework::Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{0}, "hit map on MFT disks [0,1,2,3,4]. logical-OR of each double-sided disk"};
  } dimuoncuts;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  o2::framework::HistogramRegistry fRegistry{"output", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::mt19937 engine;
  float leptonM1 = 0.f;
  float leptonM2 = 0.f;

  void init(o2::framework::InitContext& /*context*/)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());

    DefineEMEventCut();
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      DefineDielectronCut();
      leptonM1 = o2::constants::physics::MassElectron;
      leptonM2 = o2::constants::physics::MassElectron;
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      DefineDimuonCut();
      leptonM1 = o2::constants::physics::MassMuon;
      leptonM2 = o2::constants::physics::MassMuon;
    }
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

  ~DileptonProducer() {}

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

    // for pair
    fDielectronCut.SetMeeRange(dielectroncuts.cfg_min_mass, dielectroncuts.cfg_max_mass);
    fDielectronCut.SetPairPtRange(dielectroncuts.cfg_min_pair_pt, dielectroncuts.cfg_max_pair_pt);
    fDielectronCut.SetPairYRange(dielectroncuts.cfg_min_pair_y, dielectroncuts.cfg_max_pair_y);
    fDielectronCut.SetPairDCARange(dielectroncuts.cfg_min_pair_dca3d, dielectroncuts.cfg_max_pair_dca3d); // in sigma
    fDielectronCut.SetMaxMeePhiVDep([&](float phiv) { return dielectroncuts.cfg_phiv_intercept + phiv * dielectroncuts.cfg_phiv_slope; }, dielectroncuts.cfg_min_phiv, dielectroncuts.cfg_max_phiv);
    fDielectronCut.ApplyPhiV(dielectroncuts.cfg_apply_phiv);
    fDielectronCut.SetMindEtadPhi(dielectroncuts.cfg_apply_detadphi, false, dielectroncuts.cfg_min_deta, dielectroncuts.cfg_min_dphi);
    fDielectronCut.SetPairOpAng(dielectroncuts.cfg_min_opang, dielectroncuts.cfg_max_opang);
    // fDielectronCut.SetRequireDifferentSides(dielectroncuts.cfg_require_diff_sides);

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
    fDielectronCut.SetChi2TOF(0, dielectroncuts.cfg_max_chi2tof);
    // fDielectronCut.SetRelDiffPin(dielectroncuts.cfg_min_rel_diff_pin, dielectroncuts.cfg_max_rel_diff_pin);

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
    } // end of PID ML
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // for pair
    fDimuonCut.SetMassRange(dimuoncuts.cfg_min_mass, dimuoncuts.cfg_max_mass);
    fDimuonCut.SetPairPtRange(dimuoncuts.cfg_min_pair_pt, dimuoncuts.cfg_max_pair_pt);
    fDimuonCut.SetPairYRange(dimuoncuts.cfg_min_pair_y, dimuoncuts.cfg_max_pair_y);
    fDimuonCut.SetPairDCAxyRange(dimuoncuts.cfg_min_pair_dcaxy, dimuoncuts.cfg_max_pair_dcaxy);
    fDimuonCut.SetMindEtadPhi(dimuoncuts.cfg_apply_detadphi, dimuoncuts.cfg_min_deta, dimuoncuts.cfg_min_dphi);

    // for track
    fDimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, dimuoncuts.cfg_max_pt_track);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    fDimuonCut.SetTrackPhiRange(dimuoncuts.cfg_min_phi_track, dimuoncuts.cfg_max_phi_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 20);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetChi2MFT(0.f, dimuoncuts.cfg_max_chi2mft);
    // fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMaxMatchingChi2MCHMFTPtDep([&](float pt) { return (pt < dimuoncuts.cfg_border_pt_for_chi2mchmft ? dimuoncuts.cfg_max_matching_chi2_mftmch_lowPt : dimuoncuts.cfg_max_matching_chi2_mftmch_highPt); });
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
    fDimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(dimuoncuts.cfg_max_relDPt_wrt_matchedMCHMID, dimuoncuts.cfg_max_DEta_wrt_matchedMCHMID, dimuoncuts.cfg_max_DPhi_wrt_matchedMCHMID); // this is relevant for global muons
    fDimuonCut.SetMFTHitMap(dimuoncuts.requireMFTHitMap, dimuoncuts.requiredMFTDisks);
  }

  template <typename TCollision, typename TTrack1, typename TTrack2, typename TCut, typename TAllTracks>
  bool fillPairInfo(TCollision const&, TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TAllTracks const& tracks)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
        if (!cut.template IsSelectedTrack<false>(t1) || !cut.template IsSelectedTrack<false>(t2)) {
          return false;
        }
      } else { // cut-based
        if (!cut.template IsSelectedTrack<false>(t1) || !cut.template IsSelectedTrack<false>(t2)) {
          return false;
        }
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedTrack<false>(t1) || !cut.template IsSelectedTrack<false>(t2)) {
        return false;
      }

      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch(t1, cut, tracks)) {
        return false;
      }
      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch(t2, cut, tracks)) {
        return false;
      }
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.IsSelectedPair(t1, t2, d_bz, 0)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.IsSelectedPair(t1, t2)) {
        return false;
      }
    }

    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[std::make_pair(t1.globalIndex(), t2.globalIndex())];
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    // ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dca1 = 999.f, dca2 = 999.f;
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      dca1 = o2::aod::pwgem::dilepton::utils::emtrackutil::dca3DinSigma(t1);
      dca2 = o2::aod::pwgem::dilepton::utils::emtrackutil::dca3DinSigma(t2);
      if (cfgDCAType == 1) {
        dca1 = o2::aod::pwgem::dilepton::utils::emtrackutil::dcaXYinSigma(t1);
        dca2 = o2::aod::pwgem::dilepton::utils::emtrackutil::dcaXYinSigma(t2);
      } else if (cfgDCAType == 2) {
        dca1 = o2::aod::pwgem::dilepton::utils::emtrackutil::dcaZinSigma(t1);
        dca2 = o2::aod::pwgem::dilepton::utils::emtrackutil::dcaZinSigma(t2);
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      dca1 = o2::aod::pwgem::dilepton::utils::emtrackutil::fwdDcaXYinSigma(t1);
      dca2 = o2::aod::pwgem::dilepton::utils::emtrackutil::fwdDcaXYinSigma(t2);
    }

    // fill table here
    dileptonTable(eventTable.lastIndex() + 1, // lastIndex starts from -1.
                  t1.pt(), t1.eta(), t1.phi(), t1.sign(), dca1,
                  t2.pt(), t2.eta(), t2.phi(), t2.sign(), dca2,
                  weight);

    return true;
  }

  o2::framework::expressions::Filter collisionFilter_centrality = (eventcuts.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventcuts.cfgCentMax);
  o2::framework::expressions::Filter collisionFilter_numContrib = eventcuts.cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < eventcuts.cfgNumContribMax;
  o2::framework::expressions::Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  o2::framework::expressions::Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = o2::soa::Filtered<MyCollisions>;

  o2::framework::SliceCache cache;
  o2::framework::Preslice<MyElectrons> perCollision_electron = o2::aod::emprimaryelectron::emeventId;
  o2::framework::expressions::Filter trackFilter_electron = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  o2::framework::expressions::Filter ttcaFilter_electron = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  o2::framework::expressions::Filter prefilter_derived_electron = ifnode(dielectroncuts.cfg_apply_cuts_from_prefilter_derived.node() && dielectroncuts.cfg_prefilter_bits_derived.node() >= static_cast<uint16_t>(1),
                                                                         ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kMee))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kMee))) <= static_cast<uint16_t>(0), true) &&
                                                                           ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kPhiV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kPhiV))) <= static_cast<uint16_t>(0), true) &&
                                                                           ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS))) <= static_cast<uint16_t>(0), true) &&
                                                                           ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS))) <= static_cast<uint16_t>(0), true),
                                                                         o2::aod::emprimaryelectron::pfbderived >= static_cast<uint16_t>(0));

  o2::framework::expressions::Filter prefilter_electron = ifnode(dielectroncuts.cfg_apply_cuts_from_prefilter.node() && dielectroncuts.cfg_prefilter_bits.node() >= static_cast<uint16_t>(1),
                                                                 ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC))) <= static_cast<uint8_t>(0), true) &&
                                                                   ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_20MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_20MeV))) <= static_cast<uint8_t>(0), true) &&
                                                                   ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_40MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_40MeV))) <= static_cast<uint8_t>(0), true) &&
                                                                   ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_60MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_60MeV))) <= static_cast<uint8_t>(0), true) &&
                                                                   ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_80MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_80MeV))) <= static_cast<uint8_t>(0), true) &&
                                                                   ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_100MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_100MeV))) <= static_cast<uint8_t>(0), true) &&
                                                                   ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_120MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_120MeV))) <= static_cast<uint8_t>(0), true) &&
                                                                   ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_140MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_140MeV))) <= static_cast<uint8_t>(0), true),
                                                                 o2::aod::emprimaryelectron::pfb >= static_cast<uint8_t>(0));

  o2::framework::Partition<FilteredMyElectrons> positive_electrons = o2::aod::emprimaryelectron::sign > int8_t(0);
  o2::framework::Partition<FilteredMyElectrons> negative_electrons = o2::aod::emprimaryelectron::sign < int8_t(0);

  o2::framework::Preslice<MyMuons> perCollision_muon = o2::aod::emprimarymuon::emeventId;
  o2::framework::expressions::Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && o2::aod::fwdtrack::pt < dimuoncuts.cfg_max_pt_track && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track;
  o2::framework::expressions::Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);
  o2::framework::expressions::Filter prefilter_derived_muon = ifnode(dimuoncuts.cfg_apply_cuts_from_prefilter_derived.node() && dimuoncuts.cfg_prefilter_bits_derived.node() >= static_cast<uint16_t>(1),
                                                                     ifnode((dimuoncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS))) > static_cast<uint16_t>(0), (o2::aod::emprimarymuon::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS))) <= static_cast<uint16_t>(0), true) &&
                                                                       ifnode((dimuoncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS))) > static_cast<uint16_t>(0), (o2::aod::emprimarymuon::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS))) <= static_cast<uint16_t>(0), true),
                                                                     o2::aod::emprimarymuon::pfbderived >= static_cast<uint16_t>(0));
  o2::framework::Partition<FilteredMyMuons> positive_muons = o2::aod::emprimarymuon::sign > int8_t(0);
  o2::framework::Partition<FilteredMyMuons> negative_muons = o2::aod::emprimarymuon::sign < int8_t(0);

  int ndf = 0;
  template <bool isTriggerAnalysis, typename TCollisions, typename TLeptons, typename TPresilce, typename TCut, typename TAllTracks>
  void runPairing(TCollisions const& collisions, TLeptons const& posTracks, TLeptons const& negTracks, TPresilce const& perCollision, TCut const& cut, TAllTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      // float centrality = centralities[cfgCentEstimator];
      if (centralities[cfgCentEstimator] < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      float eventplanes_2_for_mix[7] = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg(), collision.ep2fv0a()};
      float ep2 = eventplanes_2_for_mix[cfgEP2Estimator_for_Mix];

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      int nuls = 0, nlspp = 0, nlsmm = 0;

      if (cfgStoreULS) {
        for (const auto& [pos, neg] : combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
          bool is_pair_ok = fillPairInfo(collision, pos, neg, cut, tracks);
          if (is_pair_ok) {
            nuls++;
          }
        }
      }
      if (cfgStoreLS) {
        for (const auto& [pos1, pos2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
          bool is_pair_ok = fillPairInfo(collision, pos1, pos2, cut, tracks);
          if (is_pair_ok) {
            nlspp++;
          }
        }
        for (const auto& [neg1, neg2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
          bool is_pair_ok = fillPairInfo(collision, neg1, neg2, cut, tracks);
          if (is_pair_ok) {
            nlsmm++;
          }
        }
      }

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        eventTable(collision.runNumber(), collision.globalBC(), collision.timestamp(), collision.posZ(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange(), centralities[cfgCentEstimator], ep2);
      }
    } // end of collision loop
  } // end of DF

  template <typename TTrack1, typename TTrack2, typename TCut, typename TAllTracks>
  bool isPairOK(TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TAllTracks const& tracks)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
        if (!cut.template IsSelectedTrack<false>(t1) || !cut.template IsSelectedTrack<false>(t2)) {
          return false;
        }
      } else { // cut-based
        if (!cut.template IsSelectedTrack<false>(t1) || !cut.template IsSelectedTrack<false>(t2)) {
          return false;
        }
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.IsSelectedTrack(t1) || !cut.IsSelectedTrack(t2)) {
        return false;
      }

      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch(t1, cut, tracks)) {
        return false;
      }
      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch(t2, cut, tracks)) {
        return false;
      }
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.IsSelectedPair(t1, t2, d_bz, 0.0)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.IsSelectedPair(t1, t2)) {
        return false;
      }
    }
    return true;
  }

  std::map<std::pair<int, int>, float> map_weight; // <posId, negId> -> float
  template <bool isTriggerAnalysis, typename TCollisions, typename TLeptons, typename TPresilce, typename TCut, typename TAllTracks>
  void fillPairWeightMap(TCollisions const& collisions, TLeptons const& posTracks, TLeptons const& negTracks, TPresilce const& perCollision, TCut const& cut, TAllTracks const& tracks)
  {
    std::vector<std::pair<int, int>> passed_pairIds;
    passed_pairIds.reserve(posTracks.size() * negTracks.size());

    for (const auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      for (const auto& [pos, neg] : combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (isPairOK(pos, neg, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(pos.globalIndex(), neg.globalIndex()));
        }
      }
      for (const auto& [pos1, pos2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (isPairOK(pos1, pos2, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(pos1.globalIndex(), pos2.globalIndex()));
        }
      }
      for (const auto& [neg1, neg2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (isPairOK(neg1, neg2, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(neg1.globalIndex(), neg2.globalIndex()));
        }
      }
    } // end of collision loop

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      for (const auto& pairId : passed_pairIds) {
        auto t1 = tracks.rawIteratorAt(std::get<0>(pairId));
        auto t2 = tracks.rawIteratorAt(std::get<1>(pairId));
        // LOGF(info, "std::get<0>(pairId) = %d, std::get<1>(pairId) = %d, t1.globalIndex() = %d, t2.globalIndex() = %d", std::get<0>(pairId), std::get<1>(pairId), t1.globalIndex(), t2.globalIndex());

        float n = 1.f; // include myself.
        for (const auto& ambId1 : t1.ambiguousElectronsIds()) {
          for (const auto& ambId2 : t2.ambiguousElectronsIds()) {
            if (std::find(passed_pairIds.begin(), passed_pairIds.end(), std::make_pair(ambId1, ambId2)) != passed_pairIds.end()) {
              n += 1.f;
            }
          }
        }
        map_weight[pairId] = 1.f / n;
      } // end of passed_pairIds loop
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      for (const auto& pairId : passed_pairIds) {
        auto t1 = tracks.rawIteratorAt(std::get<0>(pairId));
        auto t2 = tracks.rawIteratorAt(std::get<1>(pairId));

        float n = 1.f; // include myself.
        for (const auto& ambId1 : t1.ambiguousMuonsIds()) {
          for (const auto& ambId2 : t2.ambiguousMuonsIds()) {
            if (std::find(passed_pairIds.begin(), passed_pairIds.end(), std::make_pair(ambId1, ambId2)) != passed_pairIds.end()) {
              n += 1.f;
            }
          }
        }
        map_weight[pairId] = 1.f / n;
      } // end of passed_pairIds loop
    }
    passed_pairIds.clear();
    passed_pairIds.shrink_to_fit();
  }

  void processAnalysis(FilteredMyCollisions const& collisions, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<false>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons);
      }
      runPairing<false>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<false>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
      }
      runPairing<false>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
    }
    map_weight.clear();
    ndf++;
  }
  PROCESS_SWITCH(DileptonProducer, processAnalysis, "run dilepton analysis", true);

  void processNorm(o2::aod::EMEventNormInfos const& collisions)
  {
    for (const auto& collision : collisions) {
      if (collision.centFT0C() < eventcuts.cfgCentMin || eventcuts.cfgCentMax < collision.centFT0C()) {
        continue;
      }

      normTable(collision.selection_raw(), collision.rct_raw(), collision.posZint8(), collision.centFT0Muint8(), collision.centFT0Cuint8(), collision.centNTPVuint8() /*, collision.centNGlobaluint8()*/);

    } // end of collision loop
  }
  PROCESS_SWITCH(DileptonProducer, processNorm, "process normalization info", true);
};

#endif // PWGEM_DILEPTON_CORE_DILEPTONPRODUCER_H_
