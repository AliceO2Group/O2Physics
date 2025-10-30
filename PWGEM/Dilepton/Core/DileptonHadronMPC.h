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

#ifndef PWGEM_DILEPTON_CORE_DILEPTONHADRONMPC_H_
#define PWGEM_DILEPTON_CORE_DILEPTONHADRONMPC_H_

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Core/EMTrackCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMFwdTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/Dilepton/Utils/MlResponseDielectronSingleTrack.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Tools/ML/MlResponse.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"

#include "Math/Vector4D.h"
#include "TH1D.h"
#include "TString.h"

#include <algorithm>
#include <array>
#include <format>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithSWT = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMSWTriggerBits>;
using MyCollisionWithSWT = MyCollisionsWithSWT::iterator;

using MyElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit, aod::EMPrimaryElectronsPrefilterBitDerived>;
using MyElectron = MyElectrons::iterator;
using FilteredMyElectrons = soa::Filtered<MyElectrons>;
using FilteredMyElectron = FilteredMyElectrons::iterator;

using MyMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMGlobalMuonSelfIds>;
using MyMuon = MyMuons::iterator;
using FilteredMyMuons = soa::Filtered<MyMuons>;
using FilteredMyMuon = FilteredMyMuons::iterator;

using MyEMH_electron = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrack>;
using MyEMH_muon = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMFwdTrack>;
using MyEMH_track = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrack>; // for charged track

template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename TEMH, typename... Types>
struct DileptonHadronMPC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgAnalysisType{"cfgAnalysisType", static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kAzimuthalCorrelation), "kAzimuthalCorrelation:0, kCumulant:1"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth_lepton{"ndepth_lepton", 100, "depth for event mixing between lepton-lepton"};
  Configurable<int> ndepth_hadron{"ndepth_hadron", 1, "depth for event mixing between hadron-hadron"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
  Configurable<std::string> cfg_swt_name{"cfg_swt_name", "fHighTrackMult", "desired software trigger name"}; // 1 trigger name per 1 task. fHighTrackMult, fHighFt0Mult
  // Configurable<int> cfgNtracksPV08Min{"cfgNtracksPV08Min", -1, "min. multNTracksPV"};
  // Configurable<int> cfgNtracksPV08Max{"cfgNtracksPV08Max", static_cast<int>(1e+9), "max. multNTracksPV"};
  Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};
  Configurable<uint> cfgDCAType{"cfgDCAType", 0, "type of DCA for output. 0:3D, 1:XY, 2:Z, else:3D"};

  ConfigurableAxis ConfMllBins{"ConfMllBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00}, "mll bins for output histograms"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {VARIABLE_WIDTH, 0.00, 0.15, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTll bins for output histograms"};
  ConfigurableAxis ConfDCAllBins{"ConfDCAllBins", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAll bins for output histograms"};

  // ConfigurableAxis ConfMmumuBins{"ConfMmumuBins", {VARIABLE_WIDTH, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70, 5.80, 5.90, 6.00, 6.10, 6.20, 6.30, 6.40, 6.50, 6.60, 6.70, 6.80, 6.90, 7.00, 7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.80, 7.90, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.00, 11.50, 12.00}, "mmumu bins for output histograms"}; // for dimuon. one can copy bins here to hyperloop page.

  ConfigurableAxis ConfPtHadronBins{"ConfPtHadronBins", {VARIABLE_WIDTH, 0.00, 0.15, 0.2, 0.3, 0.4, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00}, "pT,h bins for output histograms"};
  ConfigurableAxis ConfYllBins{"ConfYllBins", {1, -1.f, 1.f}, "yll bins for output histograms"}; // pair rapidity
  ConfigurableAxis ConfDEtaBins{"ConfDEtaBins", {120, -6, 6}, "deta bins for output histograms"};
  Configurable<int> cfgNbinsDPhi{"cfgNbinsDPhi", 36, "nbins in dphi for output histograms"};
  Configurable<int> cfgNbinsCosNDPhi{"cfgNbinsCosNDPhi", 200, "nbins in cos(n(dphi)) for output histograms"};
  Configurable<int> cfgNmod{"cfgNmod", 2, "n-th harmonics"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
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
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } eventcuts;

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "dielectroncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pT"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pT"};
    Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -0.8, "min pair rapidity"};
    Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", +0.8, "max pair rapidity"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};
    Configurable<float> cfg_min_phiv{"cfg_min_phiv", 0.0, "min phiv (constant)"};
    Configurable<float> cfg_max_phiv{"cfg_max_phiv", 3.2, "max phiv (constant)"};
    Configurable<bool> cfg_apply_detadphi{"cfg_apply_detadphi", false, "flag to apply deta-dphi elliptic cut at PV"};
    Configurable<bool> cfg_apply_detadphiposition{"cfg_apply_detadphiposition", false, "flag to apply deta-dphi elliptic cut at certain radius"};
    Configurable<float> cfg_min_deta{"cfg_min_deta", 0.02, "min deta between 2 electrons (elliptic cut)"};
    Configurable<float> cfg_min_dphi{"cfg_min_dphi", 0.2, "min dphi between 2 electrons (elliptic cut)"};

    Configurable<bool> cfg_apply_cuts_from_prefilter{"cfg_apply_cuts_from_prefilter", false, "flag to apply prefilter set when producing derived data"};
    Configurable<uint8_t> cfg_prefilter_bits{"cfg_prefilter_bits", 0, "prefilter bits [kNone : 0, kElFromPC : 1, kElFromPi0_20MeV : 2, kElFromPi0_40MeV : 4, kElFromPi0_60MeV : 8, kElFromPi0_80MeV : 16, kElFromPi0_100MeV : 32, kElFromPi0_120MeV : 64, kElFromPi0_140MeV : 128] Please consider logical-OR among them."}; // see PairUtilities.h

    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply pair cut same as prefilter set in derived data"};
    Configurable<uint16_t> cfg_prefilter_bits_derived{"cfg_prefilter_bits_derived", 0, "prefilter bits [kNone : 0, kMee : 1, kPhiV : 2, kSplitOrMergedTrackLS : 4, kSplitOrMergedTrackULS : 8] Please consider logical-OR among them."}; // see PairUtilities.h

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
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
    Configurable<float> cfgRefR{"cfgRefR", 1.2, "reference R (in m) for extrapolation"}; // https://cds.cern.ch/record/1419204

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3, kTOFif : 4, kPIDML : 5, kTPChadrejORTOFreq_woTOFif : 6]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
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
    Configurable<bool> enableTTCA{"enableTTCA", false, "Flag to enable or disable TTCA"};

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
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pt"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pt"};
    Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -4.0, "min pair rapidity"};
    Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", -2.5, "max pair rapidity"};
    Configurable<float> cfg_min_pair_dcaxy{"cfg_min_pair_dcaxy", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dcaxy{"cfg_max_pair_dcaxy", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_detadphi{"cfg_apply_detadphi", false, "flag to apply deta-dphi elliptic cut"};
    Configurable<float> cfg_min_deta{"cfg_min_deta", 0.02, "min deta between 2 muons (elliptic cut)"};
    Configurable<float> cfg_min_dphi{"cfg_min_dphi", 0.02, "min dphi between 2 muons (elliptic cut)"};

    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 6, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 8, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+6, "max chi2"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 40, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    Configurable<bool> enableTTCA{"enableTTCA", false, "Flag to enable or disable TTCA"};
    Configurable<float> cfg_max_relDPt_wrt_matchedMCHMID{"cfg_max_relDPt_wrt_matchedMCHMID", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_DEta_wrt_matchedMCHMID{"cfg_max_DEta_wrt_matchedMCHMID", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> cfg_max_DPhi_wrt_matchedMCHMID{"cfg_max_DPhi_wrt_matchedMCHMID", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
    Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "flag to apply MFT hit map"};
    Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{0}, "hit map on MFT disks [0,1,2,3,4]. logical-OR of each double-sided disk"};
  } dimuoncuts;

  EMTrackCut fEMTrackCut;
  struct : ConfigurableGroup {
    std::string prefix = "trackcut_group";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for ref. track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 3.0, "max pT for ref. track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "min eta for ref. track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for ref. track"};
    // Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.0, "min phi for ref. track"};
    // Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for ref. track"};
    // Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.5, "max dca XY for single track in cm"};
    // Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.5, "max dca Z for single track in cm"};
    Configurable<uint16_t> cfg_track_bits{"cfg_track_bits", 5765, "required track bits"}; // default:645, loose:0, tight:778
  } trackcuts;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> cent_bin_edges;
  std::vector<float> zvtx_bin_edges;
  std::vector<float> occ_bin_edges;
  int nmod = -1;    // this is for flow analysis
  int subdet2 = -1; // this is for flow analysis
  int subdet3 = -1; // this is for flow analysis
  float leptonM1 = 0.f;
  float leptonM2 = 0.f;

  void init(InitContext& /*context*/)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);

    if (ConfVtxBins.value[0] == VARIABLE_WIDTH) {
      zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
      zvtx_bin_edges.erase(zvtx_bin_edges.begin());
      for (const auto& edge : zvtx_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: zvtx_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfVtxBins.value[0]);
      float xmin = static_cast<float>(ConfVtxBins.value[1]);
      float xmax = static_cast<float>(ConfVtxBins.value[2]);
      zvtx_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        zvtx_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: zvtx_bin_edges[%d] = %f", i, zvtx_bin_edges[i]);
      }
    }

    if (ConfCentBins.value[0] == VARIABLE_WIDTH) {
      cent_bin_edges = std::vector<float>(ConfCentBins.value.begin(), ConfCentBins.value.end());
      cent_bin_edges.erase(cent_bin_edges.begin());
      for (const auto& edge : cent_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: cent_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfCentBins.value[0]);
      float xmin = static_cast<float>(ConfCentBins.value[1]);
      float xmax = static_cast<float>(ConfCentBins.value[2]);
      cent_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        cent_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: cent_bin_edges[%d] = %f", i, cent_bin_edges[i]);
      }
    }

    LOGF(info, "cfgOccupancyEstimator = %d", cfgOccupancyEstimator.value);
    if (ConfOccupancyBins.value[0] == VARIABLE_WIDTH) {
      occ_bin_edges = std::vector<float>(ConfOccupancyBins.value.begin(), ConfOccupancyBins.value.end());
      occ_bin_edges.erase(occ_bin_edges.begin());
      for (const auto& edge : occ_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: occ_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfOccupancyBins.value[0]);
      float xmin = static_cast<float>(ConfOccupancyBins.value[1]);
      float xmax = static_cast<float>(ConfOccupancyBins.value[2]);
      occ_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        occ_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: occ_bin_edges[%d] = %f", i, occ_bin_edges[i]);
      }
    }

    emh_pos = new TEMH(ndepth_lepton);
    emh_neg = new TEMH(ndepth_lepton);
    emh_ref = new MyEMH_track(ndepth_hadron); // for reference flow

    DefineEMEventCut();
    DefineEMTrackCut();
    addhistograms();
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      DefineDielectronCut();
      leptonM1 = o2::constants::physics::MassElectron;
      leptonM2 = o2::constants::physics::MassElectron;
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      DefineDimuonCut();
      leptonM1 = o2::constants::physics::MassMuon;
      leptonM2 = o2::constants::physics::MassMuon;
    }

    if (doprocessTriggerAnalysis) {
      LOGF(info, "Trigger analysis is enabled. Desired trigger name = %s", cfg_swt_name.value);
      fRegistry.add("NormTrigger/hInspectedTVX", "inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
      fRegistry.add("NormTrigger/hScalers", "trigger counter before DS;run number;counter", kTProfile, {{80000, 520000.5, 600000.5}}, true);
      fRegistry.add("NormTrigger/hSelections", "trigger counter after DS;run number;counter", kTProfile, {{80000, 520000.5, 600000.5}}, true);
      auto hTriggerCounter = fRegistry.add<TH2>("NormTrigger/hTriggerCounter", Form("trigger counter of %s;run number;", cfg_swt_name.value.data()), kTH2D, {{80000, 520000.5, 600000.5}, {2, -0.5, 1.5}}, false);
      hTriggerCounter->GetYaxis()->SetBinLabel(1, "Analyzed Trigger");
      hTriggerCounter->GetYaxis()->SetBinLabel(2, "Analyzed TOI");
    }
  }

  template <bool isTriggerAnalysis, typename TCollision>
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
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = collision.runNumber();

    if constexpr (isTriggerAnalysis) {
      LOGF(info, "Trigger analysis is enabled. Desired trigger name = %s", cfg_swt_name.value);
      // LOGF(info, "total inspected TVX events = %d in run number %d", collision.nInspectedTVX(), collision.runNumber());
      // fRegistry.fill(HIST("Event/hNInspectedTVX"), collision.runNumber(), collision.nInspectedTVX());
    }
  }

  ~DileptonHadronMPC()
  {
    delete emh_pos;
    emh_pos = 0x0;
    delete emh_neg;
    emh_neg = 0x0;
    delete emh_ref;
    emh_ref = 0x0;

    used_trackIds_per_col.clear();
    used_trackIds_per_col.shrink_to_fit();
  }

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);

    std::string mass_axis_title = "m_{ll} (GeV/c^{2})";
    std::string pair_pt_axis_title = "p_{T,ll} (GeV/c)";
    std::string pair_dca_axis_title = "DCA_{ll} (#sigma)";
    std::string pair_rapidity_axis_title = "y_{ll}";
    std::string deta_axis_title = "#Delta#eta = #eta_{ll} - #eta_{h}";
    std::string dphi_axis_title = "#Delta#varphi = #varphi_{ll} - #varphi_{h} (rad.)";
    std::string cosndphi_axis_title = std::format("cos({0:d}(#varphi_{{ll}} - #varphi_{{h}}))", cfgNmod.value);

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      mass_axis_title = "m_{ee} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,ee} (GeV/c)";
      pair_rapidity_axis_title = "y_{ee}";
      pair_dca_axis_title = "DCA_{ee}^{3D} (#sigma)";
      if (cfgDCAType == 1) {
        pair_dca_axis_title = "DCA_{ee}^{XY} (#sigma)";
      } else if (cfgDCAType == 2) {
        pair_dca_axis_title = "DCA_{ee}^{Z} (#sigma)";
      }
      deta_axis_title = "#Delta#eta = #eta_{ee} - #eta_{h}";
      dphi_axis_title = "#Delta#varphi = #varphi_{ee} - #varphi_{h} (rad.)";
      cosndphi_axis_title = std::format("cos({0:d}(#varphi_{{ee}} - #varphi_{{h}}))", cfgNmod.value);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      mass_axis_title = "m_{#mu#mu} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,#mu#mu} (GeV/c)";
      pair_rapidity_axis_title = "y_{#mu#mu}";
      pair_dca_axis_title = "DCA_{#mu#mu}^{XY} (#sigma)";
      deta_axis_title = "#Delta#eta = #eta_{#mu#mu} - #eta_{h}";
      dphi_axis_title = "#Delta#varphi = #varphi_{#mu#mu} - #varphi_{h} (rad.)";
      cosndphi_axis_title = std::format("cos({0:d}(#varphi_{{#mu#mu}} - #varphi_{{h}}))", cfgNmod.value);
    }

    // dilepton info
    const AxisSpec axis_mass{ConfMllBins, mass_axis_title};
    const AxisSpec axis_pt{ConfPtllBins, pair_pt_axis_title};
    const AxisSpec axis_dca{ConfDCAllBins, pair_dca_axis_title};
    const AxisSpec axis_y{ConfYllBins, pair_rapidity_axis_title};

    // dilepton-hadron info
    const AxisSpec axis_deta{ConfDEtaBins, deta_axis_title};

    // hadron-hadron info
    const AxisSpec axis_deta_hh{60, -3, +3, "#Delta#eta = #eta_{h}^{ref1} - #eta_{h}^{ref2}"};

    const AxisSpec axis_pt_trg{ConfPtHadronBins, "p_{T,h} (GeV/c)"};
    const AxisSpec axis_eta_trg{40, -2, +2, "#eta_{h}"};
    const AxisSpec axis_phi_trg{36, 0, 2 * M_PI, "#varphi_{h} (rad.)"};
    fRegistry.add("Hadron/hs", "hadron", kTHnSparseD, {axis_pt_trg, axis_eta_trg, axis_phi_trg}, false);
    fRegistry.add("Hadron/hTrackBit", "track bit", kTH1D, {{65536, -0.5, 65535.5}}, false);

    fRegistry.add("Dilepton/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y}, true);
    fRegistry.addClone("Dilepton/same/uls/", "Dilepton/same/lspp/");
    fRegistry.addClone("Dilepton/same/uls/", "Dilepton/same/lsmm/");
    fRegistry.addClone("Dilepton/same/", "Dilepton/mix/");

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kAzimuthalCorrelation)) {
      const AxisSpec axis_dphi{cfgNbinsDPhi, -M_PI / 2, 3 * M_PI / 2, dphi_axis_title};
      // dilepton-hadron
      fRegistry.add("DileptonHadron/same/uls/hs", "dilepton-hadron 2PC", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y, axis_deta, axis_dphi}, true);
      fRegistry.addClone("DileptonHadron/same/uls/", "DileptonHadron/same/lspp/");
      fRegistry.addClone("DileptonHadron/same/uls/", "DileptonHadron/same/lsmm/");
      // fRegistry.addClone("DileptonHadron/same/", "DileptonHadron/mix/");

      // hadron-hadron
      const AxisSpec axis_dphi_hh{90, -M_PI / 2, 3 * M_PI / 2, "#Delta#varphi = #varphi_{h}^{ref1} - #varphi_{h}^{ref2} (rad.)"};
      fRegistry.add("HadronHadron/same/hDEtaDPhi", "hadron-hadron 2PC", kTH2D, {axis_dphi_hh, axis_deta_hh}, true);
      fRegistry.addClone("HadronHadron/same/", "HadronHadron/mix/");
      fRegistry.add("HadronHadron/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{10001, -0.5, 10000.5}}, true);
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kCumulant)) {
      const AxisSpec axis_cos_ndphi{cfgNbinsCosNDPhi, -1, +1, cosndphi_axis_title};

      // dilepton-hadron
      fRegistry.add("DileptonHadron/same/uls/hs", "dilepton-hadron 2PC", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y, axis_deta, axis_cos_ndphi}, true);
      fRegistry.addClone("DileptonHadron/same/uls/", "DileptonHadron/same/lspp/");
      fRegistry.addClone("DileptonHadron/same/uls/", "DileptonHadron/same/lsmm/");
      fRegistry.addClone("DileptonHadron/same/", "DileptonHadron/mix/");

      // hadron-hadron
      const AxisSpec axis_cosndphi_hh{cfgNbinsCosNDPhi, -1, +1, std::format("cos({0:d}(#varphi_{{h}}^{{ref1}} - #varphi_{{h}}^{{ref2}}))", cfgNmod.value)};
      fRegistry.add("HadronHadron/same/hDEtaCosNDPhi", "hadron-hadron 2PC", kTH2D, {axis_cosndphi_hh, axis_deta_hh}, true);
    }
    fRegistry.add("Dilepton/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{10001, -0.5, 10000.5}}, true);
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

    // for pair
    fDielectronCut.SetMeeRange(dielectroncuts.cfg_min_mass, dielectroncuts.cfg_max_mass);
    fDielectronCut.SetPairPtRange(dielectroncuts.cfg_min_pair_pt, dielectroncuts.cfg_max_pair_pt);
    fDielectronCut.SetPairYRange(dielectroncuts.cfg_min_pair_y, dielectroncuts.cfg_max_pair_y);
    fDielectronCut.SetPairDCARange(dielectroncuts.cfg_min_pair_dca3d, dielectroncuts.cfg_max_pair_dca3d); // in sigma
    fDielectronCut.SetMaxMeePhiVDep([&](float phiv) { return dielectroncuts.cfg_phiv_intercept + phiv * dielectroncuts.cfg_phiv_slope; }, dielectroncuts.cfg_min_phiv, dielectroncuts.cfg_max_phiv);
    fDielectronCut.ApplyPhiV(dielectroncuts.cfg_apply_phiv);
    fDielectronCut.SetMindEtadPhi(dielectroncuts.cfg_apply_detadphi, dielectroncuts.cfg_apply_detadphiposition, dielectroncuts.cfg_min_deta, dielectroncuts.cfg_min_dphi);
    fDielectronCut.SetPairOpAng(0.f, 6.3);
    fDielectronCut.SetRequireDifferentSides(false);

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, dielectroncuts.cfg_max_pt_track);
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetTrackPhiRange(dielectroncuts.cfg_min_phi_track, dielectroncuts.cfg_max_phi_track, false, false);
    fDielectronCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetMaxFracSharedClustersTPC(dielectroncuts.cfg_max_frac_shared_clusters_tpc);
    fDielectronCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITS(0.f, 16.f);
    fDielectronCut.SetTrackMaxDcaXY(dielectroncuts.cfg_max_dcaxy);
    fDielectronCut.SetTrackMaxDcaZ(dielectroncuts.cfg_max_dcaz);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);
    fDielectronCut.SetChi2TOF(0, dielectroncuts.cfg_max_chi2tof);
    fDielectronCut.SetRelDiffPin(-1e+10, +1e+10);
    fDielectronCut.IncludeITSsa(false, 0.15);

    // for eID
    fDielectronCut.SetPIDScheme(dielectroncuts.cfg_pid_scheme);
    fDielectronCut.SetTPCNsigmaElRange(dielectroncuts.cfg_min_TPCNsigmaEl, dielectroncuts.cfg_max_TPCNsigmaEl);
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
    fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
    fDimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(dimuoncuts.cfg_max_relDPt_wrt_matchedMCHMID, dimuoncuts.cfg_max_DEta_wrt_matchedMCHMID, dimuoncuts.cfg_max_DPhi_wrt_matchedMCHMID); // this is relevant for global muons
    fDimuonCut.SetMFTHitMap(dimuoncuts.requireMFTHitMap, dimuoncuts.requiredMFTDisks);
  }

  void DefineEMTrackCut()
  {
    fEMTrackCut = EMTrackCut("fEMTrackCut", "fEMTrackCut");
    fEMTrackCut.SetTrackPtRange(trackcuts.cfg_min_pt_track, trackcuts.cfg_max_pt_track);
    fEMTrackCut.SetTrackEtaRange(trackcuts.cfg_min_eta_track, trackcuts.cfg_max_eta_track);
    // fEMTrackCut.SetTrackPhiRange(trackcuts.cfg_min_phi_track, trackcuts.cfg_max_phi_track);
    // fEMTrackCut.SetTrackMaxDcaXY(trackcuts.cfg_max_dcaxy);
    // fEMTrackCut.SetTrackMaxDcaZ(trackcuts.cfg_max_dcaz);
    fEMTrackCut.SetTrackBit(trackcuts.cfg_track_bits);
  }

  template <int ev_id, typename TCollision, typename TTrack1, typename TTrack2, typename TCut, typename TAllTracks>
  bool fillDilepton(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TAllTracks const& tracks)
  {
    if constexpr (ev_id == 1) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        auto v1ambIds = t1.ambiguousElectronsIds();
        auto v2ambIds = t2.ambiguousElectronsIds();

        if ((t1.dfId() == t2.dfId()) && std::find(v2ambIds.begin(), v2ambIds.end(), t1.globalIndex()) != v2ambIds.end() && std::find(v1ambIds.begin(), v1ambIds.end(), t2.globalIndex()) != v1ambIds.end()) {
          return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
        }
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        auto v1ambIds = t1.ambiguousMuonsIds();
        auto v2ambIds = t2.ambiguousMuonsIds();

        if ((t1.dfId() == t2.dfId()) && std::find(v2ambIds.begin(), v2ambIds.end(), t1.globalIndex()) != v2ambIds.end() && std::find(v1ambIds.begin(), v1ambIds.end(), t2.globalIndex()) != v1ambIds.end()) {
          return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
        }
      }
    }

    if constexpr (ev_id == 0) {
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
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.IsSelectedPair(t1, t2, d_bz, dielectroncuts.cfgRefR)) {
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
    if (ev_id == 1) {
      weight = 1.f;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float pair_dca = 999.f;
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      pair_dca = pairDCAQuadSum(dca3DinSigma(t1), dca3DinSigma(t2));
      if (cfgDCAType == 1) {
        pair_dca = pairDCAQuadSum(dcaXYinSigma(t1), dcaXYinSigma(t2));
      } else if (cfgDCAType == 2) {
        pair_dca = pairDCAQuadSum(dcaZinSigma(t1), dcaZinSigma(t2));
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      pair_dca = pairDCAQuadSum(fwdDcaXYinSigma(t1), fwdDcaXYinSigma(t2));
    }

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Dilepton/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
      fRegistry.fill(HIST("Dilepton/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
      fRegistry.fill(HIST("Dilepton/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
    }

    // store tracks for event mixing without double counting
    if constexpr (ev_id == 0) {
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      std::vector<int> possibleIds1;
      std::vector<int> possibleIds2;

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        std::copy(t1.ambiguousElectronsIds().begin(), t1.ambiguousElectronsIds().end(), std::back_inserter(possibleIds1));
        std::copy(t2.ambiguousElectronsIds().begin(), t2.ambiguousElectronsIds().end(), std::back_inserter(possibleIds2));

        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t1.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t1.globalIndex());
          if (cfgDoMix) {
            if (t1.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMTrack(ndf, t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1, t1.cYY(), t1.cZY(), t1.cZZ()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMTrack(ndf, t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1, t1.cYY(), t1.cZY(), t1.cZZ()));
            }
          }
        }
        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t2.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t2.globalIndex());
          if (cfgDoMix) {
            if (t2.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMTrack(ndf, t2.globalIndex(), collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.dcaXY(), t2.dcaZ(), possibleIds2, t2.cYY(), t2.cZY(), t2.cZZ()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMTrack(ndf, t2.globalIndex(), collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.dcaXY(), t2.dcaZ(), possibleIds2, t2.cYY(), t2.cZY(), t2.cZZ()));
            }
          }
        }
      } else if (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        std::copy(t1.ambiguousMuonsIds().begin(), t1.ambiguousMuonsIds().end(), std::back_inserter(possibleIds1));
        std::copy(t2.ambiguousMuonsIds().begin(), t2.ambiguousMuonsIds().end(), std::back_inserter(possibleIds2));

        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t1.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t1.globalIndex());
          if (cfgDoMix) {
            if (t1.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrack(ndf, t1.globalIndex(), collision.globalIndex(), t1.fwdtrackId(), t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), possibleIds1,
                                                                        t1.cXXatDCA(), t1.cXYatDCA(), t1.cYYatDCA()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrack(ndf, t1.globalIndex(), collision.globalIndex(), t1.fwdtrackId(), t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), possibleIds1,
                                                                        t1.cXXatDCA(), t1.cXYatDCA(), t1.cYYatDCA()));
            }
          }
        }
        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t2.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t2.globalIndex());
          if (cfgDoMix) {
            if (t2.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrack(ndf, t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                        t2.cXXatDCA(), t2.cXYatDCA(), t2.cYYatDCA()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrack(ndf, t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                        t2.cXXatDCA(), t2.cXYatDCA(), t2.cYYatDCA()));
            }
          }
        }
      }
    }
    return true;
  }

  template <int ev_id, typename TTrack1, typename TTrack2, typename TCut, typename TAllTracks, typename TRefTrack>
  bool fillDileptonHadron(TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TAllTracks const& tracks, TRefTrack const& t3)
  {
    // this function must be called, if dilepton passes the cut.
    if constexpr (ev_id == 1) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        // bool is_found1 = std::find(t2.ambiguousElectronsIds.begin(), t2.ambiguousElectronsIds.end(), t1.globalIndex()) != t2.ambiguousElectronsIds.end(); // this does not work.
        // bool is_found2 = std::find(t1.ambiguousElectronsIds.begin(), t1.ambiguousElectronsIds.end(), t2.globalIndex()) != t1.ambiguousElectronsIds.end(); // this does not work.
        auto v1ambIds = t1.ambiguousElectronsIds();
        auto v2ambIds = t2.ambiguousElectronsIds();

        if ((t1.dfId() == t2.dfId()) && std::find(v2ambIds.begin(), v2ambIds.end(), t1.globalIndex()) != v2ambIds.end() && std::find(v1ambIds.begin(), v1ambIds.end(), t2.globalIndex()) != v1ambIds.end()) {
          return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
        }
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        // bool is_found1 = std::find(t2.ambiguousMuonsIds.begin(), t2.ambiguousMuonsIds.end(), t1.globalIndex()) != t2.ambiguousMuonsIds.end(); // this does not work.
        // bool is_found2 = std::find(t1.ambiguousMuonsIds.begin(), t1.ambiguousMuonsIds.end(), t2.globalIndex()) != t1.ambiguousMuonsIds.end(); // this does not work.
        auto v1ambIds = t1.ambiguousMuonsIds();
        auto v2ambIds = t2.ambiguousMuonsIds();

        if ((t1.dfId() == t2.dfId()) && std::find(v2ambIds.begin(), v2ambIds.end(), t1.globalIndex()) != v2ambIds.end() && std::find(v1ambIds.begin(), v1ambIds.end(), t2.globalIndex()) != v1ambIds.end()) {
          return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
        }
      }
    }

    if constexpr (ev_id == 0) {
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
        if (t1.trackId() == t3.trackId() || t2.trackId() == t3.trackId()) {
          return false;
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

      if (!fEMTrackCut.IsSelected(t3)) { // for charged track
        return false;
      }
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.IsSelectedPair(t1, t2, d_bz, dielectroncuts.cfgRefR)) {
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
    if (ev_id == 1) {
      weight = 1.f;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    ROOT::Math::PtEtaPhiMVector v3(t3.pt(), t3.eta(), t3.phi(), 0.139); // mass of hadron does not matter.
    float deta = v12.Eta() - v3.Eta();
    float dphi = v12.Phi() - v3.Phi();

    float pair_dca = 999.f;
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      pair_dca = pairDCAQuadSum(dca3DinSigma(t1), dca3DinSigma(t2));
      if (cfgDCAType == 1) {
        pair_dca = pairDCAQuadSum(dcaXYinSigma(t1), dcaXYinSigma(t2));
      } else if (cfgDCAType == 2) {
        pair_dca = pairDCAQuadSum(dcaZinSigma(t1), dcaZinSigma(t2));
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      pair_dca = pairDCAQuadSum(fwdDcaXYinSigma(t1), fwdDcaXYinSigma(t2));
    }

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kAzimuthalCorrelation)) {
      dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("DileptonHadron/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), deta, dphi, weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("DileptonHadron/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), deta, dphi, weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("DileptonHadron/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), deta, dphi, weight);
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kCumulant)) {
      o2::math_utils::bringTo02Pi(dphi);
      float cosndphi = std::cos(cfgNmod * dphi);
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("DileptonHadron/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), deta, cosndphi, weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("DileptonHadron/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), deta, cosndphi, weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("DileptonHadron/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), deta, cosndphi, weight);
      }
    }

    return true;
  }

  template <int ev_id, typename TRefTrack, typename TLeptons, typename TLeptonCut>
  bool fillHadronHadron(TRefTrack const& t1, TRefTrack const& t2, TLeptons const& posLeptons, TLeptons const& negLeptons, TLeptonCut const& cut)
  {
    if constexpr (ev_id == 0) {
      if (!fEMTrackCut.IsSelected(t1) || !fEMTrackCut.IsSelected(t2)) { // for charged track
        return false;
      }

      // Leptons should not be in reference track sample.
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        for (const auto& pos : posLeptons) { // leptons per collision
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut.template IsSelectedTrack<false>(pos)) {
              continue;
            }
          } else { // cut based
            if (!cut.template IsSelectedTrack<false>(pos)) {
              continue;
            }
          }
          if (t1.trackId() == pos.trackId() || t2.trackId() == pos.trackId()) {
            return false;
          }
        } // end of pos lepton loop

        for (const auto& neg : negLeptons) { // leptons per collision
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut.template IsSelectedTrack<false>(neg)) {
              continue;
            }
          } else { // cut based
            if (!cut.template IsSelectedTrack<false>(neg)) {
              continue;
            }
          }
          if (t1.trackId() == neg.trackId() || t2.trackId() == neg.trackId()) {
            return false;
          }
        } // end of neg lepton lopp

      } // end of if kDielectron
    } // end of if same event

    if constexpr (ev_id == 1) {
      if (t1.dfId() == t2.dfId() && t1.globalIndex() == t2.globalIndex()) {
        return false; // this never happens. only for protection.
      }
    }

    float weight = 1.f;
    float deta = t1.eta() - t2.eta(); // t1 is trigger, t2 is associated
    float dphi = t1.phi() - t2.phi(); // t1 is trigger, t2 is associated

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kAzimuthalCorrelation)) {
      dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
      fRegistry.fill(HIST("HadronHadron/") + HIST(event_pair_types[ev_id]) + HIST("hDEtaDPhi"), dphi, deta, weight);
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kCumulant)) {
      o2::math_utils::bringTo02Pi(dphi);
      float cosndphi = std::cos(cfgNmod * dphi);
      fRegistry.fill(HIST("HadronHadron/") + HIST(event_pair_types[ev_id]) + HIST("hDEtaCosNDPhi"), cosndphi, deta, weight);
    }
    return true;
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  // Filter collisionFilter_multiplicity = cfgNtracksPV08Min <= o2::aod::mult::multNTracksPV && o2::aod::mult::multNTracksPV < cfgNtracksPV08Max;
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  SliceCache cache;
  Preslice<MyElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Filter trackFilter_electron = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter ttcaFilter_electron = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  Filter prefilter_derived_electron = ifnode(dielectroncuts.cfg_apply_cuts_from_prefilter_derived.node() && dielectroncuts.cfg_prefilter_bits_derived.node() >= static_cast<uint16_t>(1),
                                             ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kMee))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kMee))) <= static_cast<uint16_t>(0), true) &&
                                               ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kPhiV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kPhiV))) <= static_cast<uint16_t>(0), true) &&
                                               ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS))) <= static_cast<uint16_t>(0), true) &&
                                               ifnode((dielectroncuts.cfg_prefilter_bits_derived.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfbderived & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS))) <= static_cast<uint16_t>(0), true),
                                             o2::aod::emprimaryelectron::pfbderived >= static_cast<uint16_t>(0));

  Filter prefilter_electron = ifnode(dielectroncuts.cfg_apply_cuts_from_prefilter.node() && dielectroncuts.cfg_prefilter_bits.node() >= static_cast<uint8_t>(1),
                                     ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_20MeV))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_20MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_40MeV))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_40MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_60MeV))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_60MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_80MeV))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_80MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_100MeV))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_100MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_120MeV))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_120MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_140MeV))) > static_cast<uint8_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_140MeV))) <= static_cast<uint8_t>(0), true),
                                     o2::aod::emprimaryelectron::pfb >= static_cast<uint8_t>(0));

  Partition<FilteredMyElectrons> positive_electrons = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyElectrons> negative_electrons = o2::aod::emprimaryelectron::sign < int8_t(0);

  Preslice<MyMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && o2::aod::fwdtrack::pt < dimuoncuts.cfg_max_pt_track && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track && dimuoncuts.cfg_min_phi_track < o2::aod::fwdtrack::phi && o2::aod::fwdtrack::phi < dimuoncuts.cfg_max_phi_track;
  Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);
  Partition<FilteredMyMuons> positive_muons = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<FilteredMyMuons> negative_muons = o2::aod::emprimarymuon::sign < int8_t(0);

  using RefTracks = soa::Join<aod::EMPrimaryTracks, aod::EMPrimaryTrackEMEventIds>;
  using RefTrack = RefTracks::iterator;
  Preslice<RefTracks> perCollision_track = aod::emprimarytrack::emeventId;
  Filter refTrackFilter = trackcuts.cfg_min_pt_track < 1 / nabs(o2::aod::emprimarytrack::signed1Pt) && 1 / nabs(o2::aod::emprimarytrack::signed1Pt) < trackcuts.cfg_max_pt_track && trackcuts.cfg_min_eta_track < o2::aod::emprimarytrack::eta && o2::aod::emprimarytrack::eta < trackcuts.cfg_max_eta_track;
  using FilteredRefTracks = soa::Filtered<RefTracks>;
  using FilteredRefTrack = FilteredRefTracks::iterator;

  TEMH* emh_pos = nullptr;
  TEMH* emh_neg = nullptr;
  MyEMH_track* emh_ref = nullptr; // for reference flow
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  std::vector<int> used_trackIds_per_col;
  int ndf = 0;

  template <bool isTriggerAnalysis, typename TCollisions, typename TLeptons, typename TPresilce, typename TCut, typename TAllTracks, typename TRefTracks>
  void run2PC(TCollisions const& collisions, TLeptons const& posTracks, TLeptons const& negTracks, TPresilce const& perCollision, TCut const& cut, TAllTracks const& tracks, TRefTracks const& refTracks)
  {
    for (const auto& collision : collisions) {
      initCCDB<isTriggerAnalysis>(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      float centrality = centralities[cfgCentEstimator];
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

      auto refTracks_per_coll = refTracks.sliceBy(perCollision_track, collision.globalIndex());

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      used_trackIds_per_col.reserve(posTracks_per_coll.size() + negTracks_per_coll.size());

      int nuls = 0, nlspp = 0, nlsmm = 0;
      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        bool is_pair_ok = fillDilepton<0>(collision, pos, neg, cut, tracks);
        if (is_pair_ok) {
          nuls++;
          for (const auto& refTrack : refTracks_per_coll) {
            fillDileptonHadron<0>(pos, neg, cut, tracks, refTrack);
          }
        }
      }
      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        bool is_pair_ok = fillDilepton<0>(collision, pos1, pos2, cut, tracks);
        if (is_pair_ok) {
          nlspp++;
          for (const auto& refTrack : refTracks_per_coll) {
            fillDileptonHadron<0>(pos1, pos2, cut, tracks, refTrack);
          }
        }
      }
      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        bool is_pair_ok = fillDilepton<0>(collision, neg1, neg2, cut, tracks);
        if (is_pair_ok) {
          nlsmm++;
          for (const auto& refTrack : refTracks_per_coll) {
            fillDileptonHadron<0>(neg1, neg2, cut, tracks, refTrack);
          }
        }
      }
      used_trackIds_per_col.clear();
      used_trackIds_per_col.shrink_to_fit();

      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) { // at least 1 pair exists.
        emh_ref->ReserveNTracksPerCollision(key_df_collision, refTracks_per_coll.size());
        for (const auto& track : refTracks_per_coll) {
          if (fEMTrackCut.IsSelected(track)) {
            fRegistry.fill(HIST("Hadron/hs"), track.pt(), track.eta(), track.phi());
            fRegistry.fill(HIST("Hadron/hTrackBit"), track.trackBit());

            // store ref tracks for mixed event in case of kAzimuthalCorrelation
            if (cfgDoMix && cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kAzimuthalCorrelation)) {
              emh_ref->AddTrackToEventPool(key_df_collision, EMTrack(ndf, track.globalIndex(), collision.globalIndex(), track.trackId(), track.pt(), track.eta(), track.phi(), 0.139));
            } // store ref tracks
          }
        }
        for (const auto& [ref1, ref2] : combinations(CombinationsStrictlyUpperIndexPolicy(refTracks_per_coll, refTracks_per_coll))) {
          fillHadronHadron<0>(ref1, ref2, posTracks_per_coll, negTracks_per_coll, cut);
        }
      }

      if (!cfgDoMix || !(nuls > 0 || nlspp > 0 || nlsmm > 0)) {
        continue;
      }

      // event mixing
      int zbin = lower_bound(zvtx_bin_edges.begin(), zvtx_bin_edges.end(), collision.posZ()) - zvtx_bin_edges.begin() - 1;
      if (zbin < 0) {
        zbin = 0;
      } else if (static_cast<int>(zvtx_bin_edges.size()) - 2 < zbin) {
        zbin = static_cast<int>(zvtx_bin_edges.size()) - 2;
      }

      int centbin = lower_bound(cent_bin_edges.begin(), cent_bin_edges.end(), centrality) - cent_bin_edges.begin() - 1;
      if (centbin < 0) {
        centbin = 0;
      } else if (static_cast<int>(cent_bin_edges.size()) - 2 < centbin) {
        centbin = static_cast<int>(cent_bin_edges.size()) - 2;
      }

      int epbin = 0;

      int occbin = -1;
      if (cfgOccupancyEstimator == 0) {
        occbin = lower_bound(occ_bin_edges.begin(), occ_bin_edges.end(), collision.ft0cOccupancyInTimeRange()) - occ_bin_edges.begin() - 1;
      } else if (cfgOccupancyEstimator == 1) {
        occbin = lower_bound(occ_bin_edges.begin(), occ_bin_edges.end(), collision.trackOccupancyInTimeRange()) - occ_bin_edges.begin() - 1;
      } else {
        occbin = lower_bound(occ_bin_edges.begin(), occ_bin_edges.end(), collision.ft0cOccupancyInTimeRange()) - occ_bin_edges.begin() - 1;
      }

      if (occbin < 0) {
        occbin = 0;
      } else if (static_cast<int>(occ_bin_edges.size()) - 2 < occbin) {
        occbin = static_cast<int>(occ_bin_edges.size()) - 2;
      }

      std::tuple<int, int, int, int> key_bin = std::make_tuple(zbin, centbin, epbin, occbin);

      // make a vector of selected electrons in this collision.
      auto selected_posTracks_in_this_event = emh_pos->GetTracksPerCollision(key_df_collision);
      auto selected_negTracks_in_this_event = emh_neg->GetTracksPerCollision(key_df_collision);

      auto collisionIds_in_mixing_pool = emh_pos->GetCollisionIdsFromEventPool(key_bin); // pos/neg does not matter.

      // LOGF(info, "selected_posTracks_in_this_event.size() = %d, selected_negTracks_in_this_event.size() = %d, collisionIds_in_mixing_pool.size() = %d", selected_posTracks_in_this_event.size(), selected_negTracks_in_this_event.size(), collisionIds_in_mixing_pool.size());

      // perform event mixing, only if at least 1 dilepton exists.

      for (const auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
        int mix_dfId = mix_dfId_collisionId.first;
        int mix_collisionId = mix_dfId_collisionId.second;
        if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
          continue;
        }

        auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
        uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
        fRegistry.fill(HIST("Dilepton/mix/hDiffBC"), diffBC);
        if (diffBC < ndiff_bc_mix) {
          continue;
        }

        auto posTracks_from_event_pool = emh_pos->GetTracksPerCollision(mix_dfId_collisionId);
        auto negTracks_from_event_pool = emh_neg->GetTracksPerCollision(mix_dfId_collisionId);
        // LOGF(info, "posTracks_from_event_pool.size() = %d, negTracks_from_event_pool.size() = %d", posTracks_from_event_pool.size(), negTracks_from_event_pool.size());

        for (const auto& pos : selected_posTracks_in_this_event) { // ULS mix
          for (const auto& neg : negTracks_from_event_pool) {
            fillDilepton<1>(collision, pos, neg, cut, nullptr);
          }
        }

        for (const auto& neg : selected_negTracks_in_this_event) { // ULS mix
          for (const auto& pos : posTracks_from_event_pool) {
            fillDilepton<1>(collision, neg, pos, cut, nullptr);
          }
        }

        for (const auto& pos1 : selected_posTracks_in_this_event) { // LS++ mix
          for (const auto& pos2 : posTracks_from_event_pool) {
            fillDilepton<1>(collision, pos1, pos2, cut, nullptr);
          }
        }

        for (const auto& neg1 : selected_negTracks_in_this_event) { // LS-- mix
          for (const auto& neg2 : negTracks_from_event_pool) {
            fillDilepton<1>(collision, neg1, neg2, cut, nullptr);
          }
        }
      } // end of loop over mixed event pool for lepton-lepton

      if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonHadronAnalysisType::kAzimuthalCorrelation)) {
        auto selected_refTracks_in_this_event = emh_ref->GetTracksPerCollision(key_df_collision);
        auto collisionIds_in_mixing_pool_hadron = emh_ref->GetCollisionIdsFromEventPool(key_bin);

        for (const auto& mix_dfId_collisionId : collisionIds_in_mixing_pool_hadron) {
          int mix_dfId = mix_dfId_collisionId.first;
          int mix_collisionId = mix_dfId_collisionId.second;
          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("HadronHadron/mix/hDiffBC"), diffBC);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto refTracks_from_event_pool = emh_ref->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "selected_refTracks_in_this_event.size() = %d, collisionIds_in_mixing_pool_hadron.size() = %d, refTracks_from_event_pool.size() = %d", selected_refTracks_in_this_event.size(), collisionIds_in_mixing_pool_hadron.size(), refTracks_from_event_pool.size());
          for (const auto& ref1 : selected_refTracks_in_this_event) { // ref-ref mix
            for (const auto& ref2 : refTracks_from_event_pool) {
              // LOGF(info, "ref1.pt() = %f, ref2.pt() = %f", ref1.pt(), ref2.pt());
              fillHadronHadron<1>(ref1, ref2, nullptr, nullptr, nullptr);
            }
          }
        } // end of loop over mixed event pool for hadron-hadron
      }

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
        emh_pos->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_neg->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_ref->AddCollisionIdAtLast(key_bin, key_df_collision);
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
      if (!cut.IsSelectedPair(t1, t2, d_bz, dielectroncuts.cfgRefR)) {
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
      initCCDB<isTriggerAnalysis>(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
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

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (isPairOK(pos, neg, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(pos.globalIndex(), neg.globalIndex()));
        }
      }
      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (isPairOK(pos1, pos2, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(pos1.globalIndex(), pos2.globalIndex()));
        }
      }
      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (isPairOK(neg1, neg2, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(neg1.globalIndex(), neg2.globalIndex()));
        }
      }
    } // end of collision loop

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      for (const auto& pairId : passed_pairIds) {
        auto t1 = tracks.rawIteratorAt(std::get<0>(pairId));
        auto t2 = tracks.rawIteratorAt(std::get<1>(pairId));

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

  void processAnalysis(FilteredMyCollisions const& collisions, FilteredRefTracks const& refTracks, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<false>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons);
      }
      run2PC<false>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons, refTracks);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<false>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
      }
      run2PC<false>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons, refTracks);
    }
    map_weight.clear();
    ndf++;
  }
  PROCESS_SWITCH(DileptonHadronMPC, processAnalysis, "run dilepton analysis", true);

  using FilteredMyCollisionsWithSWT = soa::Filtered<MyCollisionsWithSWT>;
  void processTriggerAnalysis(FilteredMyCollisionsWithSWT const& collisions, FilteredRefTracks const& refTracks, aod::EMSWTriggerInfos const& cefpinfos, aod::EMSWTriggerATCounters const& countersAT, aod::EMSWTriggerTOICounters const& countersTOI, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<true>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons);
      }
      run2PC<true>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons, refTracks);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<true>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
      }
      run2PC<true>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons, refTracks);
    }
    map_weight.clear();
    ndf++;

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
  PROCESS_SWITCH(DileptonHadronMPC, processTriggerAnalysis, "run dilepton analysis on triggered data", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(DileptonHadronMPC, processDummy, "Dummy function", false);
};

#endif // PWGEM_DILEPTON_CORE_DILEPTONHADRONMPC_H_
