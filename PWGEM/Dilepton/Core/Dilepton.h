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

#include <array>
#include <iterator>
#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <tuple>
#include <utility>
#include "TH1D.h"
#include "TString.h"
#include "Math/Vector4D.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "DCAFitter/DCAFitterN.h"
#include "DCAFitter/FwdDCAFitterN.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EMFwdTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithSWT = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMSWTriggerInfos>;
using MyCollisionWithSWT = MyCollisionsWithSWT::iterator;

using MyElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyElectron = MyElectrons::iterator;
using FilteredMyElectrons = soa::Filtered<MyElectrons>;
using FilteredMyElectron = FilteredMyElectrons::iterator;

using MyMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds>;
using MyMuon = MyMuons::iterator;
using FilteredMyMuons = soa::Filtered<MyMuons>;
using FilteredMyMuon = FilteredMyMuons::iterator;

using MyEMH_electron = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrackWithCov>;
using MyEMH_muon = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMFwdTrackWithCov>;

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

  Configurable<int> cfgAnalysisType{"cfgAnalysisType", static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC), "kQC:0, kUPC:1, kFlowV2:2, kFlowV3:3, kPolarization:4, kVM:5, kHFll:6"};
  Configurable<int> cfgEP2Estimator_for_Mix{"cfgEP2Estimator_for_Mix", 3, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2, NTPV:3"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 5, "difference in global BC required in mixed events"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {16, -M_PI / 2, +M_PI / 2}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
  Configurable<std::string> cfg_swt_name{"cfg_swt_name", "fHighTrackMult", "desired software trigger name"}; // 1 trigger name per 1 task. fHighTrackMult, fHighFt0Mult
  Configurable<int> cfgNtracksPV08Min{"cfgNtracksPV08Min", -1, "min. multNTracksPV"};
  Configurable<int> cfgNtracksPV08Max{"cfgNtracksPV08Max", static_cast<int>(1e+9), "max. multNTracksPV"};
  Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};
  Configurable<bool> cfgUseDCAxy{"cfgUseDCAxy", false, "flag to use DCAxy, instead of DCA3D"};

  ConfigurableAxis ConfMllBins{"ConfMllBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00}, "mll bins for output histograms"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTll bins for output histograms"};
  ConfigurableAxis ConfDCAllBins{"ConfDCAllBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAll bins for output histograms"};

  // ConfigurableAxis ConfMmumuBins{"ConfMmumuBins", {VARIABLE_WIDTH, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70, 5.80, 5.90, 6.00, 6.10, 6.20, 6.30, 6.40, 6.50, 6.60, 6.70, 6.80, 6.90, 7.00, 7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.80, 7.90, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.00, 11.50, 12.00}, "mmumu bins for output histograms"}; // for dimuon. one can copy bins here to hyperloop page.

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -1, "min. occupancy"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
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
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};

    // CCDB configuration for PID ML
    Configurable<std::string> BDTLocalPathGamma{"BDTLocalPathGamma", "pid_ml_xgboost.onnx", "Path to the local .onnx file"};
    Configurable<std::string> BDTPathCCDB{"BDTPathCCDB", "Users/d/dsekihat/pwgem/pidml/", "Path on CCDB"};
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

    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+10, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 1e+10, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
  } dimuoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  // o2::vertexing::DCAFitterN<2> fitter;
  // o2::vertexing::FwdDCAFitterN<2> fwdfitter;
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

    if (ConfVtxBins.value[0] == VARIABLE_WIDTH) {
      zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
      zvtx_bin_edges.erase(zvtx_bin_edges.begin());
      for (auto& edge : zvtx_bin_edges) {
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
      for (auto& edge : cent_bin_edges) {
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
      for (auto& edge : ep_bin_edges) {
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

    if (ConfOccupancyBins.value[0] == VARIABLE_WIDTH) {
      occ_bin_edges = std::vector<float>(ConfOccupancyBins.value.begin(), ConfOccupancyBins.value.end());
      occ_bin_edges.erase(occ_bin_edges.begin());
      for (auto& edge : occ_bin_edges) {
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
      // fitter.setPropagateToPCA(true);
      // fitter.setMaxR(5.f);
      // fitter.setMinParamChange(1e-3);
      // fitter.setMinRelChi2Change(0.9);
      // fitter.setMaxDZIni(1e9);
      // fitter.setMaxChi2(1e9);
      // fitter.setUseAbsDCA(true);
      // fitter.setWeightedFinalPCA(false);
      // fitter.setMatCorrType(matCorr);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      DefineDimuonCut();
      leptonM1 = o2::constants::physics::MassMuon;
      leptonM2 = o2::constants::physics::MassMuon;
      // fwdfitter.setPropagateToPCA(true);
      // fwdfitter.setMaxR(90.f);
      // fwdfitter.setMinParamChange(1e-3);
      // fwdfitter.setMinRelChi2Change(0.9);
      // fwdfitter.setMaxChi2(1e9);
      // fwdfitter.setUseAbsDCA(true);
      // fwdfitter.setTGeoMat(false);
    }

    fRegistry.add("Pair/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{1001, -0.5, 1000.5}}, true);
    if (doprocessTriggerAnalysis) {
      fRegistry.add("Event/hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
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
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      mRunNumber = collision.runNumber();
      // fitter.setBz(d_bz);
      // fwdfitter.setBz(d_bz);
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
    mRunNumber = collision.runNumber();
    // fitter.setBz(d_bz);
    // fwdfitter.setBz(d_bz);

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

    if constexpr (isTriggerAnalysis) {
      LOGF(info, "Trigger analysis is enabled. Desired trigger name = %s", cfg_swt_name.value);
      LOGF(info, "total inspected TVX events = %d in run number %d", collision.nInspectedTVX(), collision.runNumber());
      fRegistry.fill(HIST("Event/hNInspectedTVX"), collision.runNumber(), collision.nInspectedTVX());
    }

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

    map_mixed_eventId_to_centrality.clear();
    map_mixed_eventId_to_occupancy.clear();
    map_mixed_eventId_to_qvector.clear();
    map_mixed_eventId_to_globalBC.clear();

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();

    if (eid_bdt) {
      delete eid_bdt;
    }
    delete h2sp_resolution;
  }

  void addhistograms()
  {
    const std::string qvec_det_names[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};

    std::string mass_axis_title = "m_{ll} (GeV/c^{2})";
    std::string pair_pt_axis_title = "p_{T,ll} (GeV/c)";
    std::string pair_dca_axis_title = "DCA_{ll} (#sigma)";
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      mass_axis_title = "m_{ee} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,ee} (GeV/c)";
      pair_dca_axis_title = "DCA_{ee}^{3D} (#sigma)";
      if (cfgUseDCAxy) {
        pair_dca_axis_title = "DCA_{ee}^{XY} (#sigma)";
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      mass_axis_title = "m_{#mu#mu} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,#mu#mu} (GeV/c)";
      pair_dca_axis_title = "DCA_{#mu#mu}^{XY} (#sigma)";
    }

    // pair info
    const AxisSpec axis_mass{ConfMllBins, mass_axis_title};
    const AxisSpec axis_pt{ConfPtllBins, pair_pt_axis_title};
    const AxisSpec axis_dca{ConfDCAllBins, pair_dca_axis_title};

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC)) {
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca}, true);
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true); // phiv is only for dielectron
      }
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kUPC)) {
      const AxisSpec axis_aco{10, 0, 1.f, "#alpha = 1 - #frac{|#varphi_{l^{+}} - #varphi_{l^{-}}|}{#pi}"};
      const AxisSpec axis_asym_pt{10, 0, 1.f, "A = #frac{|p_{T,l^{+}} - p_{T,l^{-}}|}{|p_{T,l^{+}} + p_{T,l^{-}}|}"};
      const AxisSpec axis_dphi_e_ee{18, 0, M_PI, "#Delta#varphi = #varphi_{l} - #varphi_{ll} (rad.)"};
      const AxisSpec axis_cos_theta_cs{10, 0.f, 1.f, "|cos(#theta_{CS})|"};
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_aco, axis_asym_pt, axis_dphi_e_ee, axis_cos_theta_cs}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2) || cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
      if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2)) {
        nmod = 2;
      } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
        nmod = 3;
      }

      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.add("Pair/same/uls/hPrfUQ", Form("dilepton <#vec{u}_{%d,ll} #upoint #vec{Q}_{%d}^{%s}>", nmod, nmod, qvec_det_names[cfgQvecEstimator].data()), kTProfile3D, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");

      fRegistry.add("Pair/mix/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.add("Pair/mix/uls/hs_woEPmix", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.add("Pair/mix/uls/hPrfUQCosDPhi", Form("dilepton <#vec{u}_{%d,l1} #upoint #vec{Q}_{%d,ev1}^{%s} cos(%d(#varphi_{l1} - #varphi_{ll}))> + <#vec{u}_{%d,l2} #upoint #vec{Q}_{%d,ev2}^{%s} cos(%d(#varphi_{l2} - #varphi_{ll}))>", nmod, nmod, qvec_det_names[cfgQvecEstimator].data(), nmod, nmod, nmod, qvec_det_names[cfgQvecEstimator].data(), nmod), kTProfile3D, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.add("Pair/mix/uls/hPrf2UQ1UQ2CosDPhi12", Form("dilepton <2 #vec{u}_{%d,l1} #upoint #vec{Q}_{%d,ev1}^{%s} #vec{u}_{%d,l2} #upoint #vec{Q}_{%d,ev2}^{%s} cos(%d(#varphi_{l1} - #varphi_{l2}))>", nmod, nmod, qvec_det_names[cfgQvecEstimator].data(), nmod, nmod, qvec_det_names[cfgQvecEstimator].data(), nmod), kTProfile3D, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lspp/");
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lsmm/");

      // std::string detname1 = qvec_det_names[cfgQvecEstimator];
      // std::string detname2 = "";
      // std::string detname3 = "";
      // if (detname1.find("FT0") != std::string::npos) { // for dielectrons
      //   detname2 = qvec_det_names[4];
      //   detname3 = qvec_det_names[5];
      //   subdet2 = 4;
      //   subdet3 = 5;
      // } else if (detname1.find("B") != std::string::npos) { // for dimuon
      //   detname2 = qvec_det_names[1];
      //   detname3 = qvec_det_names[2];
      //   subdet2 = 1;
      //   subdet3 = 2;
      // }

      // fRegistry.add("Pair/mix/ev1/hPrf_SP12_CentFT0C", Form("Q_{%d}^{%s} #upoint Q_{%d}^{%s};centrality FT0C (%%);Q_{%d}^{%s} #upoint Q_{%d}^{%s}", nmod, detname1.data(), nmod, detname2.data(), nmod, detname1.data(), nmod, detname2.data()), kTProfile, {{110, 0, 110}}, true);
      // fRegistry.add("Pair/mix/ev1/hPrf_SP13_CentFT0C", Form("Q_{%d}^{%s} #upoint Q_{%d}^{%s};centrality FT0C (%%);Q_{%d}^{%s} #upoint Q_{%d}^{%s}", nmod, detname1.data(), nmod, detname3.data(), nmod, detname1.data(), nmod, detname3.data()), kTProfile, {{110, 0, 110}}, true);
      // fRegistry.add("Pair/mix/ev1/hPrf_SP23_CentFT0C", Form("Q_{%d}^{%s} #upoint Q_{%d}^{%s};centrality FT0C (%%);Q_{%d}^{%s} #upoint Q_{%d}^{%s}", nmod, detname2.data(), nmod, detname3.data(), nmod, detname2.data(), nmod, detname3.data()), kTProfile, {{110, 0, 110}}, true);
      // fRegistry.addClone("Pair/mix/ev1/", "Pair/mix/ev2/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kPolarization)) {
      const AxisSpec axis_cos_theta_cs{10, 0.f, 1.f, "|cos(#theta_{CS})|"};
      const AxisSpec axis_phi_cs{18, 0.f, M_PI, "|#varphi_{CS}| (rad.)"};
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_cos_theta_cs, axis_phi_cs}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kVM)) {
      std::string pair_y_axis_title = "y_{ll}";
      int nbin_y = 20;
      float min_y = -1.0;
      float max_y = +1.0;
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        pair_y_axis_title = "y_{ee}";
        nbin_y = 20;
        min_y = -1.0;
        max_y = +1.0;
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        pair_y_axis_title = "y_{#mu#mu}";
        nbin_y = 25;
        min_y = -4.5;
        max_y = -2.0;
      }
      const AxisSpec axis_y{nbin_y, min_y, max_y, pair_y_axis_title};
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_y}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kHFll)) {
      const AxisSpec axis_dphi_ee{18, 0, M_PI, "#Delta#varphi = #varphi_{e1} - #varphi_{e2} (rad.)"};
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_dphi_ee}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else { // same as kQC to avoid seg. fault
      fRegistry.add("Pair/same/uls/hs", "dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_dca}, true);
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
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetOccupancyRange(eventcuts.cfgOccupancyMin, eventcuts.cfgOccupancyMax);
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
  }

  o2::ml::OnnxModel* eid_bdt = nullptr;
  void DefineDielectronCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for pair
    fDielectronCut.SetMeeRange(dielectroncuts.cfg_min_mass, dielectroncuts.cfg_max_mass);
    fDielectronCut.SetPairPtRange(dielectroncuts.cfg_min_pair_pt, dielectroncuts.cfg_max_pair_pt);
    fDielectronCut.SetPairYRange(dielectroncuts.cfg_min_pair_y, dielectroncuts.cfg_max_pair_y);
    fDielectronCut.SetPairDCARange(dielectroncuts.cfg_min_pair_dca3d, dielectroncuts.cfg_max_pair_dca3d); // in sigma
    fDielectronCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dielectroncuts.cfg_phiv_intercept) / dielectroncuts.cfg_phiv_slope; });
    fDielectronCut.ApplyPhiV(dielectroncuts.cfg_apply_phiv);
    fDielectronCut.ApplyPrefilter(dielectroncuts.cfg_apply_pf);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, 1e+10f);
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetTrackPhiRange(dielectroncuts.cfg_min_phi_track, dielectroncuts.cfg_max_phi_track);
    fDielectronCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITS(0, 16);
    fDielectronCut.SetMaxDcaXY(dielectroncuts.cfg_max_dcaxy);
    fDielectronCut.SetMaxDcaZ(dielectroncuts.cfg_max_dcaz);

    // for eID
    fDielectronCut.SetPIDScheme(dielectroncuts.cfg_pid_scheme);
    fDielectronCut.SetTPCNsigmaElRange(dielectroncuts.cfg_min_TPCNsigmaEl, dielectroncuts.cfg_max_TPCNsigmaEl);
    fDielectronCut.SetTPCNsigmaMuRange(dielectroncuts.cfg_min_TPCNsigmaMu, dielectroncuts.cfg_max_TPCNsigmaMu);
    fDielectronCut.SetTPCNsigmaPiRange(dielectroncuts.cfg_min_TPCNsigmaPi, dielectroncuts.cfg_max_TPCNsigmaPi);
    fDielectronCut.SetTPCNsigmaKaRange(dielectroncuts.cfg_min_TPCNsigmaKa, dielectroncuts.cfg_max_TPCNsigmaKa);
    fDielectronCut.SetTPCNsigmaPrRange(dielectroncuts.cfg_min_TPCNsigmaPr, dielectroncuts.cfg_max_TPCNsigmaPr);
    fDielectronCut.SetTOFNsigmaElRange(dielectroncuts.cfg_min_TOFNsigmaEl, dielectroncuts.cfg_max_TOFNsigmaEl);

    if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      eid_bdt = new o2::ml::OnnxModel();
      if (dielectroncuts.loadModelsFromCCDB) {
        ccdbApi.init(ccdburl);
        std::map<std::string, std::string> metadata;
        bool retrieveSuccessGamma = ccdbApi.retrieveBlob(dielectroncuts.BDTPathCCDB.value, ".", metadata, dielectroncuts.timestampCCDB.value, false, dielectroncuts.BDTLocalPathGamma.value);
        if (retrieveSuccessGamma) {
          eid_bdt->initModel(dielectroncuts.BDTLocalPathGamma.value, dielectroncuts.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Gamma model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      } else {
        eid_bdt->initModel(dielectroncuts.BDTLocalPathGamma.value, dielectroncuts.enableOptimizations.value);
      }

      fDielectronCut.SetPIDModel(eid_bdt);
    } // end of PID ML
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // for pair
    fDimuonCut.SetMassRange(dimuoncuts.cfg_min_mass, dimuoncuts.cfg_max_mass);
    fDimuonCut.SetPairPtRange(dimuoncuts.cfg_min_pair_pt, dimuoncuts.cfg_max_pair_pt);
    fDimuonCut.SetPairYRange(dimuoncuts.cfg_min_pair_y, dimuoncuts.cfg_max_pair_y);
    fDimuonCut.SetPairDCAxyRange(dimuoncuts.cfg_min_pair_dcaxy, dimuoncuts.cfg_max_pair_dcaxy); // DCAxy in cm

    // for track
    fDimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, 1e10f);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_phi_track, dimuoncuts.cfg_max_phi_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 16);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
  }

  template <typename TQvectors>
  bool isGoodQvector(TQvectors const& qvectors)
  {
    bool is_good = true;
    for (auto& qn : qvectors[nmod]) {
      if (fabs(qn[0]) > 20.f || fabs(qn[1]) > 20.f) {
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

  template <int ev_id, typename TCollision, typename TTrack1, typename TTrack2, typename TCut>
  bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, TCut const& cut)
  {
    if constexpr (ev_id == 1) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        // bool is_found1 = std::find(t2.ambiguousElectronsIds.begin(), t2.ambiguousElectronsIds.end(), t1.globalIndex()) != t2.ambiguousElectronsIds.end(); //this does not work.
        // bool is_found2 = std::find(t1.ambiguousElectronsIds.begin(), t1.ambiguousElectronsIds.end(), t2.globalIndex()) != t1.ambiguousElectronsIds.end(); //this does not work.
        auto v1ambIds = t1.ambiguousElectronsIds();
        auto v2ambIds = t2.ambiguousElectronsIds();

        if ((t1.dfId() == t2.dfId()) && std::find(v2ambIds.begin(), v2ambIds.end(), t1.globalIndex()) != v2ambIds.end() && std::find(v1ambIds.begin(), v1ambIds.end(), t2.globalIndex()) != v1ambIds.end()) {
          // LOGF(info, "event id = %d: same track is found. t1.globalIndex() = %d, t1.sign() = %d, t1.pt() = %f, t1.eta() = %f, t1.phi() = %f, t2.globalIndex() = %d, t2.sign() = %d, t2.pt() = %f, t2.eta() = %f, t2.phi() = %f, deta = %f, dphi = %f (rad.)", ev_id, t1.globalIndex(), t1.sign(), t1.pt(), t1.eta(), t1.phi(), t2.globalIndex(), t2.sign(), t2.pt(), t2.eta(), t2.phi(), t1.eta() - t2.eta(), t1.phi() - t2.phi());
          return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
        }
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        // bool is_found1 = std::find(t2.ambiguousMuonsIds.begin(), t2.ambiguousMuonsIds.end(), t1.globalIndex()) != t2.ambiguousMuonsIds.end(); //this does not work.
        // bool is_found2 = std::find(t1.ambiguousMuonsIds.begin(), t1.ambiguousMuonsIds.end(), t2.globalIndex()) != t1.ambiguousMuonsIds.end(); //this does not work.
        auto v1ambIds = t1.ambiguousMuonsIds();
        auto v2ambIds = t2.ambiguousMuonsIds();

        if ((t1.dfId() == t2.dfId()) && std::find(v2ambIds.begin(), v2ambIds.end(), t1.globalIndex()) != v2ambIds.end() && std::find(v1ambIds.begin(), v1ambIds.end(), t2.globalIndex()) != v1ambIds.end()) {
          // LOGF(info, "event id = %d: same track is found. t1.globalIndex() = %d, t1.sign() = %d, t1.pt() = %f, t1.eta() = %f, t1.phi() = %f, t2.globalIndex() = %d, t2.sign() = %d, t2.pt() = %f, t2.eta() = %f, t2.phi() = %f, deta = %f, dphi = %f (rad.)", ev_id, t1.globalIndex(), t1.sign(), t1.pt(), t1.eta(), t1.phi(), t2.globalIndex(), t2.sign(), t2.pt(), t2.eta(), t2.phi(), t1.eta() - t2.eta(), t1.phi() - t2.phi());
          return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
        }
      }
    }

    if constexpr (ev_id == 0) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
          if (!cut.template IsSelectedTrack<true>(t1, collision) || !cut.template IsSelectedTrack<true>(t2, collision)) {
            return false;
          }
        } else { // cut-based
          if (!cut.template IsSelectedTrack(t1) || !cut.template IsSelectedTrack(t2)) {
            return false;
          }
        }
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        if (!cut.template IsSelectedTrack(t1) || !cut.template IsSelectedTrack(t2)) {
          return false;
        }
      }
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.template IsSelectedPair(t1, t2, d_bz)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedPair(t1, t2)) {
        return false;
      }
    }

    // float pca = 999.f, lxy = 999.f; // in unit of cm
    // o2::aod::pwgem::dilepton::utils::pairutil::isSVFound(fitter, collision, t1, t2, pca, lxy);
    // o2::aod::pwgem::dilepton::utils::pairutil::isSVFoundFwd(fwdfitter, collision, t1, t2, pca, lxy);

    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[std::make_pair(t1.globalIndex(), t2.globalIndex())];
    }
    if (ev_id == 1) {
      weight = 1.f;
    }

    // LOGF(info, "ev_id = %d, t1.sign() = %d, t2.sign() = %d, map_weight[std::make_pair(%d, %d)] = %f", ev_id, t1.sign(), t2.sign(), t1.globalIndex(), t2.globalIndex(), weight);

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dca_t1 = 999.f, dca_t2 = 999.f, pair_dca = 999.f;
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      dca_t1 = dca3DinSigma(t1);
      dca_t2 = dca3DinSigma(t2);
      pair_dca = std::sqrt((dca_t1 * dca_t1 + dca_t2 * dca_t2) / 2.);
      if (cfgUseDCAxy) {
        dca_t1 = t1.dcaXY() / std::sqrt(t1.cYY());
        dca_t2 = t2.dcaXY() / std::sqrt(t2.cYY());
        pair_dca = std::sqrt((dca_t1 * dca_t1 + dca_t2 * dca_t2) / 2.);
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      dca_t1 = fwdDcaXYinSigma(t1);
      dca_t2 = fwdDcaXYinSigma(t2);
      pair_dca = std::sqrt((dca_t1 * dca_t1 + dca_t2 * dca_t2) / 2.);
    }
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC)) {
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, weight);
        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M(), weight);
        }
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, weight);
        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M(), weight);
        }
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, weight);
        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M(), weight);
        }
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kUPC)) {
      float dphi = v1.Phi() - v2.Phi();
      o2::math_utils::bringToPMPi(dphi);
      float aco = 1.f - abs(dphi) / M_PI;
      float asym = abs(v1.Pt() - v2.Pt()) / (v1.Pt() + v2.Pt());
      float dphi_e_ee = v1.Phi() - v12.Phi();
      o2::math_utils::bringToPMPi(dphi_e_ee);

      float cos_thetaCS = 999, phiCS = 999.f;
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS<false>(t1, t2, leptonM1, leptonM2, beamE1, beamE2, beamP1, beamP2, cos_thetaCS, phiCS);

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, aco, asym, abs(dphi_e_ee), abs(cos_thetaCS), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, aco, asym, abs(dphi_e_ee), abs(cos_thetaCS), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, aco, asym, abs(dphi_e_ee), abs(cos_thetaCS), weight);
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
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, weight);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hPrfUQ"), v12.M(), v12.Pt(), pair_dca, sp, weight);
        } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, weight);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hPrfUQ"), v12.M(), v12.Pt(), pair_dca, sp, weight);
        } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, weight);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hPrfUQ"), v12.M(), v12.Pt(), pair_dca, sp, weight);
        }
      } else if constexpr (ev_id == 1) {
        if (t1.sign() * t2.sign() < 0) { // ULS
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, weight);
        } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, weight);
        } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, weight);
        }
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kPolarization)) {
      float cos_thetaCS = 999, phiCS = 999.f;
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS<false>(t1, t2, leptonM1, leptonM2, beamE1, beamE2, beamP1, beamP2, cos_thetaCS, phiCS);
      o2::math_utils::bringToPMPi(phiCS);

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, abs(cos_thetaCS), abs(phiCS), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, abs(cos_thetaCS), abs(phiCS), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, abs(cos_thetaCS), abs(phiCS), weight);
      }

    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kVM)) {
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, v12.Rapidity(), weight);
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kHFll)) {
      float dphi = v1.Phi() - v2.Phi();
      o2::math_utils::bringToPMPi(dphi);

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, abs(dphi), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, abs(dphi), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, abs(dphi), weight);
      }

    } else {                           // same as kQC to avoid seg. fault
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), pair_dca, weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), pair_dca, weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), pair_dca, weight);
      }
    }

    // store tracks for event mixing without double counting
    if constexpr (ev_id == 0) {
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());
      std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, t1.globalIndex());
      std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, t2.globalIndex());

      std::vector<int> possibleIds1;
      std::vector<int> possibleIds2;

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        std::copy(t1.ambiguousElectronsIds().begin(), t1.ambiguousElectronsIds().end(), std::back_inserter(possibleIds1));
        std::copy(t2.ambiguousElectronsIds().begin(), t2.ambiguousElectronsIds().end(), std::back_inserter(possibleIds2));

        if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id1) == used_trackIds.end()) {
          used_trackIds.emplace_back(pair_tmp_id1);
          if (cfgDoMix) {
            if (t1.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMTrackWithCov(ndf, t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1,
                                                                            t1.x(), t1.y(), t1.z(), t1.alpha(), t1.snp(), t1.tgl(), t1.cYY(), t1.cZY(), t1.cZZ(),
                                                                            t1.cSnpY(), t1.cSnpZ(), t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(), t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMTrackWithCov(ndf, t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), leptonM1, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1,
                                                                            t1.x(), t1.y(), t1.z(), t1.alpha(), t1.snp(), t1.tgl(), t1.cYY(), t1.cZY(), t1.cZZ(),
                                                                            t1.cSnpY(), t1.cSnpZ(), t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(), t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()));
            }
          }
        }
        if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id2) == used_trackIds.end()) {
          used_trackIds.emplace_back(pair_tmp_id2);
          if (cfgDoMix) {
            if (t2.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMTrackWithCov(ndf, t2.globalIndex(), collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.dcaXY(), t2.dcaZ(), possibleIds2,
                                                                            t2.x(), t2.y(), t2.z(), t2.alpha(), t2.snp(), t2.tgl(), t2.cYY(), t2.cZY(), t2.cZZ(),
                                                                            t2.cSnpY(), t2.cSnpZ(), t2.cSnpSnp(), t2.cTglY(), t2.cTglZ(), t2.cTglSnp(), t2.cTglTgl(), t2.c1PtY(), t2.c1PtZ(), t2.c1PtSnp(), t2.c1PtTgl(), t2.c1Pt21Pt2()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMTrackWithCov(ndf, t2.globalIndex(), collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), leptonM2, t2.sign(), t2.dcaXY(), t2.dcaZ(), possibleIds2,
                                                                            t2.x(), t2.y(), t2.z(), t2.alpha(), t2.snp(), t2.tgl(), t2.cYY(), t2.cZY(), t2.cZZ(),
                                                                            t2.cSnpY(), t2.cSnpZ(), t2.cSnpSnp(), t2.cTglY(), t2.cTglZ(), t2.cTglSnp(), t2.cTglTgl(), t2.c1PtY(), t2.c1PtZ(), t2.c1PtSnp(), t2.c1PtTgl(), t2.c1Pt21Pt2()));
            }
          }
        }
      } else if (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        std::copy(t1.ambiguousMuonsIds().begin(), t1.ambiguousMuonsIds().end(), std::back_inserter(possibleIds1));
        std::copy(t2.ambiguousMuonsIds().begin(), t2.ambiguousMuonsIds().end(), std::back_inserter(possibleIds2));

        if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id1) == used_trackIds.end()) {
          used_trackIds.emplace_back(pair_tmp_id1);
          if (cfgDoMix) {
            if (t1.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(ndf, t1.globalIndex(), collision.globalIndex(), t1.fwdtrackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), possibleIds1,
                                                                               t1.x(), t1.y(), t1.z(), t1.tgl(), t1.cXX(), t1.cXY(), t1.cYY(),
                                                                               t1.cPhiX(), t1.cPhiY(), t1.cPhiPhi(), t1.cTglX(), t1.cTglY(), t1.cTglPhi(), t1.cTglTgl(), t1.c1PtX(), t1.c1PtY(), t1.c1PtPhi(), t1.c1PtTgl(), t1.c1Pt21Pt2(), t1.chi2()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(ndf, t1.globalIndex(), collision.globalIndex(), t1.fwdtrackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), possibleIds1,
                                                                               t1.x(), t1.y(), t1.z(), t1.tgl(), t1.cXX(), t1.cXY(), t1.cYY(),
                                                                               t1.cPhiX(), t1.cPhiY(), t1.cPhiPhi(), t1.cTglX(), t1.cTglY(), t1.cTglPhi(), t1.cTglTgl(), t1.c1PtX(), t1.c1PtY(), t1.c1PtPhi(), t1.c1PtTgl(), t1.c1Pt21Pt2(), t1.chi2()));
            }
          }
        }
        if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id2) == used_trackIds.end()) {
          used_trackIds.emplace_back(pair_tmp_id2);
          if (cfgDoMix) {
            if (t2.sign() > 0) {
              emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(ndf, t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                               t2.x(), t2.y(), t2.z(), t2.tgl(), t2.cXX(), t2.cXY(), t2.cYY(),
                                                                               t2.cPhiX(), t2.cPhiY(), t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(), t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2(), t2.chi2()));
            } else {
              emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(ndf, t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                               t2.x(), t2.y(), t2.z(), t2.tgl(), t2.cXX(), t2.cXY(), t2.cYY(),
                                                                               t2.cPhiX(), t2.cPhiY(), t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(), t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2(), t2.chi2()));
            }
          }
        }
      }
    }
    return true;
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_multiplicity = cfgNtracksPV08Min <= o2::aod::mult::multNTracksPV && o2::aod::mult::multNTracksPV < cfgNtracksPV08Max;
  Filter collisionFilter_occupancy = eventcuts.cfgOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgOccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  SliceCache cache;
  Preslice<MyElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Filter trackFilter_electron = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && dielectroncuts.cfg_min_phi_track < o2::aod::track::phi && o2::aod::track::phi < dielectroncuts.cfg_max_phi_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter_electron = (dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < dielectroncuts.cfg_min_TPCNsigmaPi || dielectroncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi);
  Filter ttcaFilter_electron = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  Partition<FilteredMyElectrons> positive_electrons = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyElectrons> negative_electrons = o2::aod::emprimaryelectron::sign < int8_t(0);

  Preslice<MyMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track && dimuoncuts.cfg_min_phi_track < o2::aod::fwdtrack::phi && o2::aod::fwdtrack::phi < dimuoncuts.cfg_max_phi_track;
  Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);
  Partition<FilteredMyMuons> positive_muons = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<FilteredMyMuons> negative_muons = o2::aod::emprimarymuon::sign < int8_t(0);

  TEMH* emh_pos = nullptr;
  TEMH* emh_neg = nullptr;
  std::map<std::pair<int, int>, std::vector<std::vector<std::array<float, 2>>>> map_mixed_eventId_to_qvector;
  std::map<std::pair<int, int>, float> map_mixed_eventId_to_centrality;
  std::map<std::pair<int, int>, int> map_mixed_eventId_to_occupancy;
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  std::vector<std::pair<int, int>> used_trackIds;
  int ndf = 0;

  template <bool isTriggerAnalysis, typename TCollisions, typename TLeptons, typename TPresilce, typename TCut>
  void runPairing(TCollisions const& collisions, TLeptons const& posTracks, TLeptons const& negTracks, TPresilce const& perCollision, TCut const& cut)
  {
    for (auto& collision : collisions) {
      initCCDB<isTriggerAnalysis>(collision);
      const float centralities[4] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()};
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

      int nuls = 0, nlspp = 0, nlsmm = 0;
      for (auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        bool is_pair_ok = fillPairInfo<0>(collision, pos, neg, cut);
        if (is_pair_ok) {
          nuls++;
        }
      }
      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        bool is_pair_ok = fillPairInfo<0>(collision, pos1, pos2, cut);
        if (is_pair_ok) {
          nlspp++;
        }
      }
      for (auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        bool is_pair_ok = fillPairInfo<0>(collision, neg1, neg2, cut);
        if (is_pair_ok) {
          nlsmm++;
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

      int epbin = lower_bound(ep_bin_edges.begin(), ep_bin_edges.end(), ep2) - ep_bin_edges.begin() - 1;
      if (epbin < 0) {
        epbin = 0;
      } else if (static_cast<int>(ep_bin_edges.size()) - 2 < epbin) {
        epbin = static_cast<int>(ep_bin_edges.size()) - 2;
      }

      int occbin = lower_bound(occ_bin_edges.begin(), occ_bin_edges.end(), collision.trackOccupancyInTimeRange()) - occ_bin_edges.begin() - 1;
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

      for (auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
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

        for (auto& pos : selected_posTracks_in_this_event) { // ULS mix
          for (auto& neg : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos, neg, cut);
          }
        }

        for (auto& neg : selected_negTracks_in_this_event) { // ULS mix
          for (auto& pos : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, neg, pos, cut);
          }
        }

        for (auto& pos1 : selected_posTracks_in_this_event) { // LS++ mix
          for (auto& pos2 : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos1, pos2, cut);
          }
        }

        for (auto& neg1 : selected_negTracks_in_this_event) { // LS-- mix
          for (auto& neg2 : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, neg1, neg2, cut);
          }
        }
      } // end of loop over mixed event pool

      // run mixed event loop for flow measurement. Don't divide mixed-event categories by event planes, if you do flow measurement.
      if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2) || cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
        for (int epbin_tmp = 0; epbin_tmp < static_cast<int>(ep_bin_edges.size()) - 1; epbin_tmp++) {
          std::tuple<int, int, int, int> key_bin = std::make_tuple(zbin, centbin, epbin_tmp, occbin);
          auto collisionIds_in_mixing_pool = emh_pos->GetCollisionIdsFromEventPool(key_bin); // pos/neg does not matter.
          // LOGF(info, "collisionIds_in_mixing_pool.size() = %d", collisionIds_in_mixing_pool.size());

          for (auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
            int mix_dfId = mix_dfId_collisionId.first;
            int mix_collisionId = mix_dfId_collisionId.second;
            if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
              continue;
            }

            auto centrality_mix = map_mixed_eventId_to_centrality[mix_dfId_collisionId];
            auto occupancy_mix = map_mixed_eventId_to_occupancy[mix_dfId_collisionId];
            auto qvectors_mix = map_mixed_eventId_to_qvector[mix_dfId_collisionId];
            auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
            uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
            if (diffBC < ndiff_bc_mix) {
              continue;
            }

            auto posTracks_from_event_pool = emh_pos->GetTracksPerCollision(mix_dfId_collisionId);
            auto negTracks_from_event_pool = emh_neg->GetTracksPerCollision(mix_dfId_collisionId);
            // LOGF(info, "Do event mixing: current event (%d, %d) | event pool (%d, %d), npos = %d , nneg = %d", ndf, collision.globalIndex(), mix_dfId, mix_collisionId, posTracks_from_event_pool.size(), negTracks_from_event_pool.size());

            for (auto& pos : selected_posTracks_in_this_event) { // ULS mix
              for (auto& neg : negTracks_from_event_pool) {
                fillMixedPairInfoForFlow(pos, neg, cut, qvectors, qvectors_mix, collision.centFT0C(), centrality_mix, collision.trackOccupancyInTimeRange(), occupancy_mix);
              }
            }

            for (auto& neg : selected_negTracks_in_this_event) { // ULS mix
              for (auto& pos : posTracks_from_event_pool) {
                fillMixedPairInfoForFlow(neg, pos, cut, qvectors, qvectors_mix, collision.centFT0C(), centrality_mix, collision.trackOccupancyInTimeRange(), occupancy_mix);
              }
            }

            for (auto& pos1 : selected_posTracks_in_this_event) { // LS++ mix
              for (auto& pos2 : posTracks_from_event_pool) {
                fillMixedPairInfoForFlow(pos1, pos2, cut, qvectors, qvectors_mix, collision.centFT0C(), centrality_mix, collision.trackOccupancyInTimeRange(), occupancy_mix);
              }
            }

            for (auto& neg1 : selected_negTracks_in_this_event) { // LS-- mix
              for (auto& neg2 : negTracks_from_event_pool) {
                fillMixedPairInfoForFlow(neg1, neg2, cut, qvectors, qvectors_mix, collision.centFT0C(), centrality_mix, collision.trackOccupancyInTimeRange(), occupancy_mix);
              }
            }
          } // end of loop over mixed event pool

        } // end of ep bin loop

      } // end of mixing loop for flow

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        if (nmod > 0) {
          map_mixed_eventId_to_qvector[key_df_collision] = qvectors;
          map_mixed_eventId_to_centrality[key_df_collision] = collision.centFT0C();
          map_mixed_eventId_to_occupancy[key_df_collision] = collision.trackOccupancyInTimeRange();
        }
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
        emh_pos->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_neg->AddCollisionIdAtLast(key_bin, key_df_collision);
      }

    } // end of collision loop

    ndf++;
  } // end of DF

  template <typename TTrack1, typename TTrack2, typename TCut, typename TQvectors, typename TMixedQvectors>
  bool fillMixedPairInfoForFlow(TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TQvectors const& qvectors, TMixedQvectors const& qvectors_mix, const float centrality, const float centrality_mix, const int occupancy, const int occupancy_mix)
  {
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

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.template IsSelectedPair(t1, t2, d_bz)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedPair(t1, t2)) {
        return false;
      }
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dca_t1 = 999.f, dca_t2 = 999.f, pair_dca = 999.f;
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      dca_t1 = dca3DinSigma(t1);
      dca_t2 = dca3DinSigma(t2);
      pair_dca = std::sqrt((dca_t1 * dca_t1 + dca_t2 * dca_t2) / 2.);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      dca_t1 = fwdDcaXYinSigma(t1);
      dca_t2 = fwdDcaXYinSigma(t2);
      pair_dca = std::sqrt((dca_t1 * dca_t1 + dca_t2 * dca_t2) / 2.);
    }

    float sp1 = RecoDecay::dotProd(std::array<float, 2>{static_cast<float>(std::cos(nmod * v1.Phi())), static_cast<float>(std::sin(nmod * v1.Phi()))}, qvectors[nmod][cfgQvecEstimator]) / getSPresolution(centrality, occupancy);
    float sp2 = RecoDecay::dotProd(std::array<float, 2>{static_cast<float>(std::cos(nmod * v2.Phi())), static_cast<float>(std::sin(nmod * v2.Phi()))}, qvectors_mix[nmod][cfgQvecEstimator]) / getSPresolution(centrality_mix, occupancy_mix);
    float cos_dphi1 = std::cos(nmod * (v1.Phi() - v12.Phi()));
    float cos_dphi2 = std::cos(nmod * (v2.Phi() - v12.Phi()));
    float cos_dphi12 = std::cos(nmod * (v1.Phi() - v2.Phi()));
    const float weight = 1.f;

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Pair/mix/uls/hs_woEPmix"), v12.M(), v12.Pt(), pair_dca, weight);
      fRegistry.fill(HIST("Pair/mix/uls/hPrfUQCosDPhi"), v12.M(), v12.Pt(), pair_dca, sp1 * cos_dphi1 + sp2 * cos_dphi2, weight);
      fRegistry.fill(HIST("Pair/mix/uls/hPrf2UQ1UQ2CosDPhi12"), v12.M(), v12.Pt(), pair_dca, 2.f * sp1 * sp2 * cos_dphi12, weight);
    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
      fRegistry.fill(HIST("Pair/mix/lspp/hs_woEPmix"), v12.M(), v12.Pt(), pair_dca, weight);
      fRegistry.fill(HIST("Pair/mix/lspp/hPrfUQCosDPhi"), v12.M(), v12.Pt(), pair_dca, sp1 * cos_dphi1 + sp2 * cos_dphi2, weight);
      fRegistry.fill(HIST("Pair/mix/lspp/hPrf2UQ1UQ2CosDPhi12"), v12.M(), v12.Pt(), pair_dca, 2.f * sp1 * sp2 * cos_dphi12, weight);
    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
      fRegistry.fill(HIST("Pair/mix/lsmm/hs_woEPmix"), v12.M(), v12.Pt(), pair_dca, weight);
      fRegistry.fill(HIST("Pair/mix/lsmm/hPrfUQCosDPhi"), v12.M(), v12.Pt(), pair_dca, sp1 * cos_dphi1 + sp2 * cos_dphi2, weight);
      fRegistry.fill(HIST("Pair/mix/lsmm/hPrf2UQ1UQ2CosDPhi12"), v12.M(), v12.Pt(), pair_dca, 2.f * sp1 * sp2 * cos_dphi12, weight);
    }
    return true;
  }

  template <typename TCollision, typename TTrack1, typename TTrack2, typename TCut>
  bool isPairOK(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, TCut const& cut)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
        if (!cut.template IsSelectedTrack<true>(t1, collision) || !cut.template IsSelectedTrack<true>(t2, collision)) {
          return false;
        }
      } else { // cut-based
        if (!cut.template IsSelectedTrack(t1) || !cut.template IsSelectedTrack(t2)) {
          return false;
        }
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedTrack(t1) || !cut.template IsSelectedTrack(t2)) {
        return false;
      }
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.template IsSelectedPair(t1, t2, d_bz)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedPair(t1, t2)) {
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

    for (auto& collision : collisions) {
      initCCDB<isTriggerAnalysis>(collision);
      const float centralities[4] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()};
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
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 4th harmonics
      };

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      if (nmod > 0 && !isGoodQvector(qvectors)) {
        continue;
      }

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      for (auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (isPairOK(collision, pos, neg, cut)) {
          passed_pairIds.emplace_back(std::make_pair(pos.globalIndex(), neg.globalIndex()));
        }
      }
      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (isPairOK(collision, pos1, pos2, cut)) {
          passed_pairIds.emplace_back(std::make_pair(pos1.globalIndex(), pos2.globalIndex()));
        }
      }
      for (auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (isPairOK(collision, neg1, neg2, cut)) {
          passed_pairIds.emplace_back(std::make_pair(neg1.globalIndex(), neg2.globalIndex()));
        }
      }
    } // end of collision loop

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      for (auto& pairId : passed_pairIds) {
        auto t1 = tracks.rawIteratorAt(std::get<0>(pairId));
        auto t2 = tracks.rawIteratorAt(std::get<1>(pairId));
        // LOGF(info, "std::get<0>(pairId) = %d, std::get<1>(pairId) = %d, t1.globalIndex() = %d, t2.globalIndex() = %d", std::get<0>(pairId), std::get<1>(pairId), t1.globalIndex(), t2.globalIndex());

        float n = 1.f; // include myself.
        for (auto& ambId1 : t1.ambiguousElectronsIds()) {
          for (auto& ambId2 : t2.ambiguousElectronsIds()) {
            if (std::find(passed_pairIds.begin(), passed_pairIds.end(), std::make_pair(ambId1, ambId2)) != passed_pairIds.end()) {
              n += 1.f;
            }
          }
        }
        map_weight[pairId] = 1.f / n;
      } // end of passed_pairIds loop
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      for (auto& pairId : passed_pairIds) {
        auto t1 = tracks.rawIteratorAt(std::get<0>(pairId));
        auto t2 = tracks.rawIteratorAt(std::get<1>(pairId));

        float n = 1.f; // include myself.
        for (auto& ambId1 : t1.ambiguousMuonsIds()) {
          for (auto& ambId2 : t2.ambiguousMuonsIds()) {
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
      runPairing<false>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<false>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
      }
      runPairing<false>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut);
    }
    map_weight.clear();
  }
  PROCESS_SWITCH(Dilepton, processAnalysis, "run dilepton analysis", true);

  using FilteredMyCollisionsWithSWT = soa::Filtered<MyCollisionsWithSWT>;
  void processTriggerAnalysis(FilteredMyCollisionsWithSWT const& collisions, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<true>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, electrons);
      }
      runPairing<true>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap<true>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, muons);
      }
      runPairing<true>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut);
    }
    map_weight.clear();
  }
  PROCESS_SWITCH(Dilepton, processTriggerAnalysis, "run dilepton analysis on triggered data", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(Dilepton, processDummy, "Dummy function", false);
};

#endif // PWGEM_DILEPTON_CORE_DILEPTON_H_
