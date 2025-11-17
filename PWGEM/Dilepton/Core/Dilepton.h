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

#ifndef PWGEM_DILEPTON_CORE_DILEPTON_H_
#define PWGEM_DILEPTON_CORE_DILEPTON_H_

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
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

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithSWT = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMSWTriggerBits>;
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

template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename TEMH, typename... Types>
struct Dilepton {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<bool> cfgApplySPresolution{"cfgApplySPresolution", false, "flag to apply resolution correction for flow analysis"};
  Configurable<std::string> spresoPath{"spresoPath", "Users/d/dsekihat/PWGEM/dilepton/Qvector/resolution/LHC23zzh/pass3/test", "Path to SP resolution file"};
  Configurable<std::string> spresoHistName{"spresoHistName", "h1_R2_FT0M_BPos_BNeg", "histogram name of SP resolution file"};

  Configurable<int> cfgAnalysisType{"cfgAnalysisType", static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC), "kQC:0, kUPC:1, kFlowV2:2, kFlowV3:3, kPolarization:4, kHFll:5"};
  Configurable<int> cfgEP2Estimator_for_Mix{"cfgEP2Estimator_for_Mix", 3, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {16, -M_PI / 2, +M_PI / 2}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
  Configurable<std::string> cfg_swt_name{"cfg_swt_name", "fHighTrackMult", "desired software trigger name"}; // 1 trigger name per 1 task. fHighTrackMult, fHighFt0Mult
  Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
  Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
  Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};
  Configurable<uint> cfgDCAType{"cfgDCAType", 0, "type of DCA for output. 0:3D, 1:XY, 2:Z, else:3D"};
  Configurable<bool> cfgUseSignedDCA{"cfgUseSignedDCA", false, "flag to use signs in the DCA calculation"};
  Configurable<int> cfgPolarizationFrame{"cfgPolarizationFrame", 0, "frame of polarization. 0:CS, 1:HX, else:FATAL"};

  ConfigurableAxis ConfMllBins{"ConfMllBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00}, "mll bins for output histograms"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTll bins for output histograms"};
  ConfigurableAxis ConfDCAllBins{"ConfDCAllBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAll bins for output histograms"};
  ConfigurableAxis ConfYllBins{"ConYllBins", {1, -1.f, 1.f}, "yll bins for output histograms"}; // pair rapidity

  // ConfigurableAxis ConfMmumuBins{"ConfMmumuBins", {VARIABLE_WIDTH, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70, 5.80, 5.90, 6.00, 6.10, 6.20, 6.30, 6.40, 6.50, 6.60, 6.70, 6.80, 6.90, 7.00, 7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.80, 7.90, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.00, 11.50, 12.00}, "mmumu bins for output histograms"}; // for dimuon. one can copy bins here to hyperloop page.

  ConfigurableAxis ConfSPBins{"ConfSPBins", {200, -5, 5}, "SP bins for flow analysis"};
  ConfigurableAxis ConfPolarizationCosThetaBins{"ConfPolarizationCosThetaBins", {20, -1.f, 1.f}, "cos(theta) bins for polarization analysis"};
  ConfigurableAxis ConfPolarizationPhiBins{"ConfPolarizationPhiBins", {1, -M_PI, M_PI}, "phi bins for polarization analysis"};
  ConfigurableAxis ConfPolarizationQuadMomBins{"ConfPolarizationQuadMomBins", {15, -0.5, 1}, "quadrupole moment bins for polarization analysis"}; // quardrupole moment <(3 x cos^2(theta) -1)/2>

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
    Configurable<float> cfg_min_deta{"cfg_min_deta", 0.02, "min deta between 2 electrons (elliptic cut)"};
    Configurable<float> cfg_min_dphi{"cfg_min_dphi", 0.2, "min dphi between 2 electrons (elliptic cut)"};
    Configurable<float> cfg_min_opang{"cfg_min_opang", 0.0, "min opening angle"};
    Configurable<float> cfg_max_opang{"cfg_max_opang", 6.4, "max opening angle"};
    Configurable<bool> cfg_require_diff_sides{"cfg_require_diff_sides", false, "flag to require 2 tracks are from different sides."};

    Configurable<bool> cfg_apply_cuts_from_prefilter{"cfg_apply_cuts_from_prefilter", false, "flag to apply prefilter set when producing derived data"};
    Configurable<uint16_t> cfg_prefilter_bits{"cfg_prefilter_bits", 0, "prefilter bits [kNone : 0, kElFromPC : 1, kElFromPi0_20MeV : 2, kElFromPi0_40MeV : 4, kElFromPi0_60MeV : 8, kElFromPi0_80MeV : 16, kElFromPi0_100MeV : 32, kElFromPi0_120MeV : 64, kElFromPi0_140MeV : 128] Please consider logical-OR among them."}; // see PairUtilities.h

    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply pair cut same as prefilter set in derived data"};
    Configurable<uint16_t> cfg_prefilter_bits_derived{"cfg_prefilter_bits_derived", 0, "prefilter bits [kNone : 0, kMee : 1, kPhiV : 2, kSplitOrMergedTrackLS : 4, kSplitOrMergedTrackULS : 8] Please consider logical-OR among them."}; // see PairUtilities.h

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
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

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3, kTOFif : 4, kPIDML : 5, kTPChadrejORTOFreq_woTOFif : 6]"};
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
  int mRunNumber;
  float d_bz;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> cent_bin_edges;
  std::vector<float> zvtx_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;
  int nmod = -1;    // this is for flow analysis
  int subdet2 = -1; // this is for flow analysis
  int subdet3 = -1; // this is for flow analysis
  float leptonM1 = 0.f;
  float leptonM2 = 0.f;

  float beamM1 = o2::constants::physics::MassProton; // mass of beam
  float beamM2 = o2::constants::physics::MassProton; // mass of beam
  float beamE1 = 0.f;                                // beam energy
  float beamE2 = 0.f;                                // beam energy
  float beamP1 = 0.f;                                // beam momentum
  float beamP2 = 0.f;                                // beam momentum
  TH2D* h2sp_resolution = nullptr;

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

    if (ConfEPBins.value[0] == VARIABLE_WIDTH) {
      ep_bin_edges = std::vector<float>(ConfEPBins.value.begin(), ConfEPBins.value.end());
      ep_bin_edges.erase(ep_bin_edges.begin());
      for (const auto& edge : ep_bin_edges) {
        LOGF(info, "VARIABLE_WIDTH: ep_bin_edges = %f", edge);
      }
    } else {
      int nbins = static_cast<int>(ConfEPBins.value[0]);
      float xmin = static_cast<float>(ConfEPBins.value[1]);
      float xmax = static_cast<float>(ConfEPBins.value[2]);
      ep_bin_edges.resize(nbins + 1);
      for (int i = 0; i < nbins + 1; i++) {
        ep_bin_edges[i] = (xmax - xmin) / (nbins)*i + xmin;
        LOGF(info, "FIXED_WIDTH: ep_bin_edges[%d] = %f", i, ep_bin_edges[i]);
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

    emh_pos = new TEMH(ndepth);
    emh_neg = new TEMH(ndepth);

    DefineEMEventCut();
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

    fRegistry.add("Pair/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{10001, -0.5, 10000.5}}, true);

    if (doprocessNorm) {
      fRegistry.addClone("Event/before/hCollisionCounter", "Event/norm/hCollisionCounter");
    }
    if (doprocessTriggerAnalysis) {
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

    auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", collision.timestamp());
    int beamZ1 = grplhcif->getBeamZ(o2::constants::lhc::BeamC);
    int beamZ2 = grplhcif->getBeamZ(o2::constants::lhc::BeamA);
    int beamA1 = grplhcif->getBeamA(o2::constants::lhc::BeamC);
    int beamA2 = grplhcif->getBeamA(o2::constants::lhc::BeamA);
    beamE1 = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamC);
    beamE2 = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamA);
    beamM1 = o2::constants::physics::MassProton * beamA1;
    beamM2 = o2::constants::physics::MassProton * beamA2;
    beamP1 = std::sqrt(std::pow(beamE1, 2) - std::pow(beamM1, 2));
    beamP2 = std::sqrt(std::pow(beamE2, 2) - std::pow(beamM2, 2));
    LOGF(info, "beamZ1 = %d, beamZ2 = %d, beamA1 = %d, beamA2 = %d, beamE1 = %f (GeV), beamE2 = %f (GeV), beamM1 = %f (GeV), beamM2 = %f (GeV), beamP1 = %f (GeV), beamP2 = %f (GeV)", beamZ1, beamZ2, beamA1, beamA2, beamE1, beamE2, beamM1, beamM2, beamP1, beamP2);

    if (cfgApplySPresolution) {
      auto list = ccdb->getForTimeStamp<TList>(spresoPath, collision.timestamp());
      h2sp_resolution = reinterpret_cast<TH2D*>(list->FindObject(spresoHistName.value.data()));
      LOGF(info, "h2sp_resolution.GetBinContent(40, 1) = %f", h2sp_resolution->GetBinContent(40, 1));
    }
  }

  ~Dilepton()
  {
    delete emh_pos;
    emh_pos = 0x0;
    delete emh_neg;
    emh_neg = 0x0;

    used_trackIds_per_col.clear();
    used_trackIds_per_col.shrink_to_fit();
    map_mixed_eventId_to_globalBC.clear();

    delete h2sp_resolution;
  }

  void addhistograms()
  {
    const std::string qvec_det_names[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};

    std::string mass_axis_title = "m_{ll} (GeV/c^{2})";
    std::string pair_pt_axis_title = "p_{T,ll} (GeV/c)";
    std::string pair_dca_axis_title = "DCA_{ll} (#sigma)";
    std::string pair_y_axis_title = "y_{ll}";
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      mass_axis_title = "m_{ee} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,ee} (GeV/c)";
      pair_y_axis_title = "y_{ee}";
      pair_dca_axis_title = "DCA_{ee}^{3D} (#sigma)";
      if (cfgDCAType == 1) {
        pair_dca_axis_title = "DCA_{ee}^{XY} (#sigma)";
      } else if (cfgDCAType == 2) {
        pair_dca_axis_title = "DCA_{ee}^{Z} (#sigma)";
      }
      if (cfgUseSignedDCA) {
        pair_dca_axis_title = "Signed " + pair_dca_axis_title;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      mass_axis_title = "m_{#mu#mu} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,#mu#mu} (GeV/c)";
      pair_dca_axis_title = "DCA_{#mu#mu}^{XY} (#sigma)";
      pair_y_axis_title = "y_{#mu#mu}";
    }

    // pair info
    const AxisSpec axis_mass{ConfMllBins, mass_axis_title};
    const AxisSpec axis_pt{ConfPtllBins, pair_pt_axis_title};
    const AxisSpec axis_dca{ConfDCAllBins, pair_dca_axis_title};
    const AxisSpec axis_y{ConfYllBins, pair_y_axis_title};

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC)) {
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y}, true);
      fRegistry.add("Pair/same/uls/hDeltaEtaDeltaPhi", "#Delta#eta-#Delta#varphi between 2 tracks;#Delta#varphi (rad.);#Delta#eta;", kTH2D, {{180, -M_PI, M_PI}, {400, -2, +2}}, true);

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 1.0f}}, true); // phiv is only for dielectron
        fRegistry.add("Pair/same/uls/hMvsOpAng", "m_{ee} vs. angle between 2 tracks;#omega (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 1.0f}}, true);
        fRegistry.add("Pair/same/uls/hDCA1vsDCA2", "DCA of leg1 vs. DCA of leg2;DCA1(#sigma);DCA2 (#sigma)", kTH2D, {{100, 0, 10.0}, {100, 0, 10}}, true);
      }
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kUPC)) {
      const AxisSpec axis_aco{10, 0, 1.f, "#alpha = 1 - #frac{|#varphi_{l^{+}} - #varphi_{l^{-}}|}{#pi}"};
      const AxisSpec axis_asym_pt{10, 0, 1.f, "A = #frac{|p_{T,l^{+}} - p_{T,l^{-}}|}{|p_{T,l^{+}} + p_{T,l^{-}}|}"};
      const AxisSpec axis_dphi_e_ee{18, 0, M_PI, "#Delta#varphi = #varphi_{l} - #varphi_{ll} (rad.)"};

      std::string frameName = "CS";
      if (cfgPolarizationFrame == 0) {
        frameName = "CS";
      } else if (cfgPolarizationFrame == 1) {
        frameName = "HX";
      } else {
        LOG(fatal) << "set 0 or 1 to cfgPolarizationFrame!";
      }
      const AxisSpec axis_cos_theta{ConfPolarizationCosThetaBins, Form("cos(#theta^{%s})", frameName.data())};
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y, axis_aco, axis_asym_pt, axis_dphi_e_ee, axis_cos_theta}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2) || cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
      if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2)) {
        nmod = 2;
      } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
        nmod = 3;
      }

      const AxisSpec axis_sp{ConfSPBins, Form("#vec{u}_{%d,ll} #upoint #vec{Q}_{%d}^{%s}", nmod, nmod, qvec_det_names[cfgQvecEstimator].data())};

      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y, axis_sp}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");

      fRegistry.add("Pair/mix/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y}, true);
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lspp/");
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lsmm/");

    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kPolarization)) {
      std::string frameName = "CS";
      if (cfgPolarizationFrame == 0) {
        frameName = "CS";
      } else if (cfgPolarizationFrame == 1) {
        frameName = "HX";
      } else {
        LOG(fatal) << "set 0 or 1 to cfgPolarizationFrame!";
      }

      const AxisSpec axis_cos_theta{ConfPolarizationCosThetaBins, Form("cos(#theta^{%s})", frameName.data())};
      const AxisSpec axis_phi{ConfPolarizationPhiBins, Form("#varphi^{%s} (rad.)", frameName.data())};
      const AxisSpec axis_quadmom{ConfPolarizationQuadMomBins, Form("#frac{3 cos^{2}(#theta^{%s}) -1}{2}", frameName.data())};
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y, axis_cos_theta, axis_phi, axis_quadmom}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kHFll)) {
      const AxisSpec axis_dphi_ee{36, -M_PI / 2., 3. / 2. * M_PI, "#Delta#varphi = #varphi_{l1} - #varphi_{l2} (rad.)"}; // for kHFll
      const AxisSpec axis_deta_ee{40, -2., 2., "#Delta#eta = #eta_{l1} - #eta_{l2}"};
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y, axis_dphi_ee, axis_deta_ee}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else { // same as kQC to avoid seg. fault
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    }

    // event info
    if (nmod == 2) {
      o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<2>(&fRegistry);
    } else if (nmod == 3) {
      o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<3>(&fRegistry);
    } else {
      o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);
    }
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix", Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEP2Estimator_for_Mix].data()), kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry.add("Event/after/hEP2_CentFT0C_forMix", Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEP2Estimator_for_Mix].data()), kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
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
    fDielectronCut.SetMindEtadPhi(dielectroncuts.cfg_apply_detadphi, false, dielectroncuts.cfg_min_deta, dielectroncuts.cfg_min_dphi);
    fDielectronCut.SetPairOpAng(dielectroncuts.cfg_min_opang, dielectroncuts.cfg_max_opang);
    fDielectronCut.SetRequireDifferentSides(dielectroncuts.cfg_require_diff_sides);

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
    fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
    fDimuonCut.SetMaxdPtdEtadPhiwrtMCHMID(dimuoncuts.cfg_max_relDPt_wrt_matchedMCHMID, dimuoncuts.cfg_max_DEta_wrt_matchedMCHMID, dimuoncuts.cfg_max_DPhi_wrt_matchedMCHMID); // this is relevant for global muons
    fDimuonCut.SetMFTHitMap(dimuoncuts.requireMFTHitMap, dimuoncuts.requiredMFTDisks);
  }

  template <typename TQvectors>
  bool isGoodQvector(TQvectors const& qvectors)
  {
    bool is_good = true;
    for (const auto& qn : qvectors[nmod]) {
      if (std::fabs(qn[0]) > 20.f || std::fabs(qn[1]) > 20.f) {
        is_good = false;
        break;
      }
    }
    return is_good;
  }

  float getSPresolution(const float centrality, const int occupancy)
  {
    if (h2sp_resolution == nullptr) {
      return 1.f;
    } else {
      int binId_cen = h2sp_resolution->GetXaxis()->FindBin(centrality);
      int binId_occ = h2sp_resolution->GetYaxis()->FindBin(occupancy);

      if (centrality < h2sp_resolution->GetXaxis()->GetXmin()) {
        binId_cen = 1;
      }
      if (h2sp_resolution->GetXaxis()->GetXmax() < centrality) {
        binId_cen = h2sp_resolution->GetXaxis()->GetNbins();
      }

      if (occupancy < h2sp_resolution->GetYaxis()->GetXmin()) {
        binId_occ = 1;
      }
      if (h2sp_resolution->GetYaxis()->GetXmax() < occupancy) {
        binId_occ = h2sp_resolution->GetYaxis()->GetNbins();
      }
      return h2sp_resolution->GetBinContent(binId_cen, binId_occ);
    }
  }

  template <int ev_id, typename TCollision, typename TTrack1, typename TTrack2, typename TCut, typename TAllTracks>
  bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TAllTracks const& tracks)
  {
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
      if (!cut.IsSelectedPair(t1, t2, d_bz, 0)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.IsSelectedPair(t1, t2)) {
        return false;
      }
    }

    float weight = 1.f;
    if constexpr (ev_id == 0) {
      if (cfgApplyWeightTTCA) {
        weight = map_weight[std::make_pair(t1.globalIndex(), t2.globalIndex())];
      }
    }
    // LOGF(info, "ev_id = %d, t1.sign() = %d, t2.sign() = %d, map_weight[std::make_pair(%d, %d)] = %f", ev_id, t1.sign(), t2.sign(), t1.globalIndex(), t2.globalIndex(), weight);

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float pair_dca = 999.f;
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (cfgUseSignedDCA) {
        pair_dca = pairDCASignQuadSum(dca3DinSigma(t1), dca3DinSigma(t2), t1.sign(), t2.sign());
      } else {
        pair_dca = pairDCAQuadSum(dca3DinSigma(t1), dca3DinSigma(t2));
      }
      if (cfgDCAType == 1) {
        if (cfgUseSignedDCA) {
          pair_dca = pairDCASignQuadSum(dcaXYinSigma(t1), dcaXYinSigma(t2), t1.sign(), t2.sign());
        } else {
          pair_dca = pairDCAQuadSum(dcaXYinSigma(t1), dcaXYinSigma(t2));
        }
      } else if (cfgDCAType == 2) {
        if (cfgUseSignedDCA) {
          pair_dca = pairDCASignQuadSum(dcaZinSigma(t1), dcaZinSigma(t2), t1.sign(), t2.sign());
        } else {
          pair_dca = pairDCAQuadSum(dcaZinSigma(t1), dcaZinSigma(t2));
        }
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      pair_dca = pairDCAQuadSum(fwdDcaXYinSigma(t1), fwdDcaXYinSigma(t2));
    }

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC)) {
      float deta = t1.sign() * v1.Pt() > t2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
      float dphi = t1.sign() * v1.Pt() > t2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
      o2::math_utils::bringToPMPi(dphi);

      float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
      float opAng = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hDeltaEtaDeltaPhi"), dphi, deta, weight);
        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M(), weight);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsOpAng"), opAng, v12.M(), weight);
          if (cfgDCAType == 1) {
            fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hDCA1vsDCA2"), dcaXYinSigma(t1), dcaXYinSigma(t2), weight);
          } else if (cfgDCAType == 2) {
            fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hDCA1vsDCA2"), dcaZinSigma(t1), dcaZinSigma(t2), weight);
          }
        }
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hDeltaEtaDeltaPhi"), dphi, deta, weight);
        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M(), weight);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsOpAng"), opAng, v12.M(), weight);
          if (cfgDCAType == 1) {
            fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hDCA1vsDCA2"), dcaXYinSigma(t1), dcaXYinSigma(t2), weight);
          } else if (cfgDCAType == 2) {
            fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hDCA1vsDCA2"), dcaZinSigma(t1), dcaZinSigma(t2), weight);
          }
        }
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hDeltaEtaDeltaPhi"), dphi, deta, weight);
        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M(), weight);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsOpAng"), opAng, v12.M(), weight);
          if (cfgDCAType == 1) {
            fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hDCA1vsDCA2"), dcaXYinSigma(t1), dcaXYinSigma(t2), weight);
          } else if (cfgDCAType == 2) {
            fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hDCA1vsDCA2"), dcaZinSigma(t1), dcaZinSigma(t2), weight);
          }
        }
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kUPC)) {
      float dphi = v1.Phi() - v2.Phi();
      o2::math_utils::bringToPMPi(dphi);
      float aco = 1.f - std::fabs(dphi) / M_PI;
      float asym = std::fabs(v1.Pt() - v2.Pt()) / (v1.Pt() + v2.Pt());
      float dphi_e_ee = v1.Phi() - v12.Phi();
      o2::math_utils::bringToPMPi(dphi_e_ee);

      float cos_thetaPol = 999, phiPol = 999.f;
      if (cfgPolarizationFrame == 0) {
        o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(std::array<float, 4>{t1.px(), t1.py(), t1.pz(), leptonM1}, std::array<float, 4>{t2.px(), t2.py(), t2.pz(), leptonM2}, beamE1, beamE2, beamP1, beamP2, t1.sign(), cos_thetaPol, phiPol);
      } else if (cfgPolarizationFrame == 1) {
        o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(std::array<float, 4>{t1.px(), t1.py(), t1.pz(), leptonM1}, std::array<float, 4>{t2.px(), t2.py(), t2.pz(), leptonM2}, beamE1, beamE2, beamP1, beamP2, t1.sign(), cos_thetaPol, phiPol);
      }
      o2::math_utils::bringToPMPi(phiPol);

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), aco, asym, std::fabs(dphi_e_ee), cos_thetaPol, weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), aco, asym, std::fabs(dphi_e_ee), cos_thetaPol, weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), aco, asym, std::fabs(dphi_e_ee), cos_thetaPol, weight);
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2) || cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2btot = {collision.q2xbtot(), collision.q2ybtot()};
      std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
      std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};
      std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
      std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
      std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
      std::array<float, 2> q3btot = {collision.q3xbtot(), collision.q3ybtot()};
      std::array<float, 2> q3bpos = {collision.q3xbpos(), collision.q3ybpos()};
      std::array<float, 2> q3bneg = {collision.q3xbneg(), collision.q3ybneg()};

      std::vector<std::vector<std::array<float, 2>>> qvectors = {
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 0th harmonics
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 1st harmonics
        {q2ft0m, q2ft0a, q2ft0c, q2btot, q2bpos, q2bneg},                                                 // 2nd harmonics
        {q3ft0m, q3ft0a, q3ft0c, q3btot, q3bpos, q3bneg},                                                 // 3rd harmonics
      };

      if constexpr (ev_id == 0) {
        // LOGF(info, "collision.centFT0C() = %f, collision.trackOccupancyInTimeRange() = %d, getSPresolution = %f", collision.centFT0C(), collision.trackOccupancyInTimeRange(), getSPresolution(collision.centFT0C(), collision.trackOccupancyInTimeRange()));

        float sp = RecoDecay::dotProd(std::array<float, 2>{static_cast<float>(std::cos(nmod * v12.Phi())), static_cast<float>(std::sin(nmod * v12.Phi()))}, qvectors[nmod][cfgQvecEstimator]) / getSPresolution(collision.centFT0C(), collision.trackOccupancyInTimeRange());
        if (t1.sign() * t2.sign() < 0) { // ULS
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), sp, weight);
        } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), sp, weight);
        } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), sp, weight);
        }
      } else if constexpr (ev_id == 1) {
        if (t1.sign() * t2.sign() < 0) { // ULS
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
        } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
        } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
        }
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kPolarization)) {
      float cos_thetaPol = 999, phiPol = 999.f;
      if (cfgPolarizationFrame == 0) {
        o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(std::array<float, 4>{t1.px(), t1.py(), t1.pz(), leptonM1}, std::array<float, 4>{t2.px(), t2.py(), t2.pz(), leptonM2}, beamE1, beamE2, beamP1, beamP2, t1.sign(), cos_thetaPol, phiPol);
      } else if (cfgPolarizationFrame == 1) {
        o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(std::array<float, 4>{t1.px(), t1.py(), t1.pz(), leptonM1}, std::array<float, 4>{t2.px(), t2.py(), t2.pz(), leptonM2}, beamE1, beamE2, beamP1, beamP2, t1.sign(), cos_thetaPol, phiPol);
      }
      o2::math_utils::bringToPMPi(phiPol);
      float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), cos_thetaPol, phiPol, quadmom, weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), cos_thetaPol, phiPol, quadmom, weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), cos_thetaPol, phiPol, quadmom, weight);
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kHFll)) {
      float dphi = v1.Phi() - v2.Phi();
      dphi = RecoDecay::constrainAngle(dphi, -o2::constants::math::PIHalf);
      float deta = v1.Eta() - v2.Eta();

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), dphi, deta, weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), dphi, deta, weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), dphi, deta, weight);
      }

    } else {                           // same as kQC to avoid seg. fault
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
      }
    }

    // store tracks for event mixing without double counting
    if constexpr (ev_id == 0) {
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t1.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t1.globalIndex());
          if (cfgDoMix) {
            if (t1.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMTrack(t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.dcaXY(), t1.dcaZ(), t1.cYY(), t1.cZY(), t1.cZZ()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMTrack(t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.dcaXY(), t1.dcaZ(), t1.cYY(), t1.cZY(), t1.cZZ()));
            }
          }
        }
        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t2.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t2.globalIndex());
          if (cfgDoMix) {
            if (t2.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMTrack(t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.dcaXY(), t2.dcaZ(), t2.cYY(), t2.cZY(), t2.cZZ()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMTrack(t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.dcaXY(), t2.dcaZ(), t2.cYY(), t2.cZY(), t2.cZZ()));
            }
          }
        }
      } else if (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t1.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t1.globalIndex());
          if (cfgDoMix) {
            if (t1.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrack(t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), t1.cXXatDCA(), t1.cXYatDCA(), t1.cYYatDCA()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrack(t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), t1.cXXatDCA(), t1.cXYatDCA(), t1.cYYatDCA()));
            }
          }
        }
        if (std::find(used_trackIds_per_col.begin(), used_trackIds_per_col.end(), t2.globalIndex()) == used_trackIds_per_col.end()) {
          used_trackIds_per_col.emplace_back(t2.globalIndex());
          if (cfgDoMix) {
            if (t2.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrack(t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), t2.cXXatDCA(), t2.cXYatDCA(), t2.cYYatDCA()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrack(t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), t2.cXXatDCA(), t2.cXYatDCA(), t2.cYYatDCA()));
            }
          }
        }
      }
    }
    return true;
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_numContrib = cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < cfgNumContribMax;
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

  Filter prefilter_electron = ifnode(dielectroncuts.cfg_apply_cuts_from_prefilter.node() && dielectroncuts.cfg_prefilter_bits.node() >= static_cast<uint16_t>(1),
                                     ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_20MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_20MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_40MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_40MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_60MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_60MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_80MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_80MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_100MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_100MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_120MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_120MeV))) <= static_cast<uint8_t>(0), true) &&
                                       ifnode((dielectroncuts.cfg_prefilter_bits.node() & static_cast<uint16_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_140MeV))) > static_cast<uint16_t>(0), (o2::aod::emprimaryelectron::pfb & static_cast<uint8_t>(1 << int(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_140MeV))) <= static_cast<uint8_t>(0), true),
                                     o2::aod::emprimaryelectron::pfb >= static_cast<uint8_t>(0));

  Partition<FilteredMyElectrons> positive_electrons = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyElectrons> negative_electrons = o2::aod::emprimaryelectron::sign < int8_t(0);

  Preslice<MyMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && o2::aod::fwdtrack::pt < dimuoncuts.cfg_max_pt_track && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track;
  Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);
  Partition<FilteredMyMuons> positive_muons = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<FilteredMyMuons> negative_muons = o2::aod::emprimarymuon::sign < int8_t(0);

  TEMH* emh_pos = nullptr;
  TEMH* emh_neg = nullptr;
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  std::vector<int> used_trackIds_per_col;
  int ndf = 0;

  template <bool isTriggerAnalysis, typename TCollisions, typename TLeptons, typename TPresilce, typename TCut, typename TAllTracks>
  void runPairing(TCollisions const& collisions, TLeptons const& posTracks, TLeptons const& negTracks, TPresilce const& perCollision, TCut const& cut, TAllTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
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

      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2btot = {collision.q2xbtot(), collision.q2ybtot()};
      std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
      std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};
      std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
      std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
      std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
      std::array<float, 2> q3btot = {collision.q3xbtot(), collision.q3ybtot()};
      std::array<float, 2> q3bpos = {collision.q3xbpos(), collision.q3ybpos()};
      std::array<float, 2> q3bneg = {collision.q3xbneg(), collision.q3ybneg()};
      const float eventplanes_2_for_mix[6] = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      float ep2 = eventplanes_2_for_mix[cfgEP2Estimator_for_Mix];

      std::vector<std::vector<std::array<float, 2>>> qvectors = {
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 0th harmonics
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 1st harmonics
        {q2ft0m, q2ft0a, q2ft0c, q2btot, q2bpos, q2bneg},                                                 // 2nd harmonics
        {q3ft0m, q3ft0a, q3ft0c, q3btot, q3bpos, q3bneg},                                                 // 3rd harmonics
      };

      if (nmod == 2) {
        o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, 2>(&fRegistry, collision);
      } else if (nmod == 3) {
        o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, 3>(&fRegistry, collision);
      } else {
        o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision);
      }
      if (nmod < 0 || isGoodQvector(qvectors)) {
        fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev - 1); // is qvector calibarated
      }
      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      if (nmod > 0 && !isGoodQvector(qvectors)) {
        continue;
      }

      if (nmod == 2) {
        o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, 2>(&fRegistry, collision);
      } else if (nmod == 3) {
        o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, 3>(&fRegistry, collision);
      } else {
        o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision);
      }
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev - 1); // is qvector calibarated

      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      // LOGF(info, "collision.globalIndex() = %d , collision.posZ() = %f , collision.numContrib() = %d, centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", collision.globalIndex(), collision.posZ(), collision.numContrib(), centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      used_trackIds_per_col.reserve(posTracks_per_coll.size() + negTracks_per_coll.size());
      int nuls = 0, nlspp = 0, nlsmm = 0;
      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        bool is_pair_ok = fillPairInfo<0>(collision, pos, neg, cut, tracks);
        if (is_pair_ok) {
          nuls++;
        }
      }
      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        bool is_pair_ok = fillPairInfo<0>(collision, pos1, pos2, cut, tracks);
        if (is_pair_ok) {
          nlspp++;
        }
      }
      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        bool is_pair_ok = fillPairInfo<0>(collision, neg1, neg2, cut, tracks);
        if (is_pair_ok) {
          nlsmm++;
        }
      }
      used_trackIds_per_col.clear();
      used_trackIds_per_col.shrink_to_fit();

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

      int epbin = lower_bound(ep_bin_edges.begin(), ep_bin_edges.end(), ep2) - ep_bin_edges.begin() - 1;
      if (epbin < 0) {
        epbin = 0;
      } else if (static_cast<int>(ep_bin_edges.size()) - 2 < epbin) {
        epbin = static_cast<int>(ep_bin_edges.size()) - 2;
      }

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

      // LOGF(info, "collision.globalIndex() = %d, collision.posZ() = %f, centrality = %f, ep2 = %f, collision.trackOccupancyInTimeRange() = %d, zbin = %d, centbin = %d, epbin = %d, occbin = %d", collision.globalIndex(), collision.posZ(), centrality, ep2, collision.trackOccupancyInTimeRange(), zbin, centbin, epbin, occbin);

      std::tuple<int, int, int, int> key_bin = std::make_tuple(zbin, centbin, epbin, occbin);
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      // make a vector of selected photons in this collision.
      auto selected_posTracks_in_this_event = emh_pos->GetTracksPerCollision(key_df_collision);
      auto selected_negTracks_in_this_event = emh_neg->GetTracksPerCollision(key_df_collision);
      // LOGF(info, "N selected tracks in current event (%d, %d), zvtx = %f, centrality = %f , npos = %d , nneg = %d, nuls = %d , nlspp = %d, nlsmm = %d", ndf, collision.globalIndex(), collision.posZ(), centralities[cfgCentEstimator], selected_posTracks_in_this_event.size(), selected_negTracks_in_this_event.size(), nuls, nlspp, nlsmm);

      auto collisionIds_in_mixing_pool = emh_pos->GetCollisionIdsFromEventPool(key_bin); // pos/neg does not matter.
      // LOGF(info, "collisionIds_in_mixing_pool.size() = %d", collisionIds_in_mixing_pool.size());

      for (const auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
        int mix_dfId = mix_dfId_collisionId.first;
        int mix_collisionId = mix_dfId_collisionId.second;
        if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
          continue;
        }

        auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
        uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
        fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
        if (diffBC < ndiff_bc_mix) {
          continue;
        }

        auto posTracks_from_event_pool = emh_pos->GetTracksPerCollision(mix_dfId_collisionId);
        auto negTracks_from_event_pool = emh_neg->GetTracksPerCollision(mix_dfId_collisionId);
        // LOGF(info, "Do event mixing: current event (%d, %d) | event pool (%d, %d), npos = %d , nneg = %d", ndf, collision.globalIndex(), mix_dfId, mix_collisionId, posTracks_from_event_pool.size(), negTracks_from_event_pool.size());

        for (const auto& pos : selected_posTracks_in_this_event) { // ULS mix
          for (const auto& neg : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos, neg, cut, nullptr);
          }
        }

        for (const auto& neg : selected_negTracks_in_this_event) { // ULS mix
          for (const auto& pos : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, neg, pos, cut, nullptr);
          }
        }

        for (const auto& pos1 : selected_posTracks_in_this_event) { // LS++ mix
          for (const auto& pos2 : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos1, pos2, cut, nullptr);
          }
        }

        for (const auto& neg1 : selected_negTracks_in_this_event) { // LS-- mix
          for (const auto& neg2 : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, neg1, neg2, cut, nullptr);
          }
        }
      } // end of loop over mixed event pool

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
        emh_pos->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_neg->AddCollisionIdAtLast(key_bin, key_df_collision);
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
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (!collision.swtalias_bit(o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value))) {
          continue;
        }
      }

      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2btot = {collision.q2xbtot(), collision.q2ybtot()};
      std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
      std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};
      std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
      std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
      std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
      std::array<float, 2> q3btot = {collision.q3xbtot(), collision.q3ybtot()};
      std::array<float, 2> q3bpos = {collision.q3xbpos(), collision.q3ybpos()};
      std::array<float, 2> q3bneg = {collision.q3xbneg(), collision.q3ybneg()};

      std::vector<std::vector<std::array<float, 2>>> qvectors = {
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 0th harmonics
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 1st harmonics
        {q2ft0m, q2ft0a, q2ft0c, q2btot, q2bpos, q2bneg},                                                 // 2nd harmonics
        {q3ft0m, q3ft0a, q3ft0c, q3btot, q3bpos, q3bneg},                                                 // 3rd harmonics
      };

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      if (nmod > 0 && !isGoodQvector(qvectors)) {
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
  PROCESS_SWITCH(Dilepton, processAnalysis, "run dilepton analysis", true);

  using FilteredMyCollisionsWithSWT = soa::Filtered<MyCollisionsWithSWT>;
  void processTriggerAnalysis(FilteredMyCollisionsWithSWT const& collisions, aod::EMSWTriggerInfos const& cefpinfos, aod::EMSWTriggerATCounters const& countersAT, aod::EMSWTriggerTOICounters const& countersTOI, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<true>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons);
      }
      runPairing<true>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<true>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
      }
      runPairing<true>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
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
  PROCESS_SWITCH(Dilepton, processTriggerAnalysis, "run dilepton analysis on triggered data", false);

  Filter collisionFilter_centrality_norm = cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax;
  using FilteredNormInfos = soa::Filtered<aod::EMEventNormInfos>;

  void processNorm(FilteredNormInfos const& collisions)
  {
    for (const auto& collision : collisions) {
      if (collision.centFT0C() < cfgCentMin || cfgCentMax < collision.centFT0C()) {
        continue;
      }

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
  PROCESS_SWITCH(Dilepton, processNorm, "process normalization info", false);

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
  PROCESS_SWITCH(Dilepton, processBC, "process BC counter", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(Dilepton, processDummy, "Dummy function", false);
};

#endif // PWGEM_DILEPTON_CORE_DILEPTON_H_
