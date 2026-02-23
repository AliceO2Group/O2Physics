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
// This code runs loop over leptons in MC
//    Please write to: daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_CORE_DILEPTONMC_H_
#define PWGEM_DILEPTON_CORE_DILEPTONMC_H_

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
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

#include "Math/Vector4D.h"
#include "TString.h"

#include <algorithm>
#include <array>
#include <format>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::MostProbableEMEventIdsInMC>;
using MyMCCollision = MyMCCollisions::iterator;

using MyMCElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit, aod::EMPrimaryElectronsPrefilterBitDerived, aod::EMPrimaryElectronMCLabels>;
using MyMCElectron = MyMCElectrons::iterator;
using FilteredMyMCElectrons = soa::Filtered<MyMCElectrons>;
using FilteredMyMCElectron = FilteredMyMCElectrons::iterator;

using MyMCMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMGlobalMuonSelfIds, aod::EMPrimaryMuonMCLabels, aod::EMMFTMCLabels>;
using MyMCMuon = MyMCMuons::iterator;
using FilteredMyMCMuons = soa::Filtered<MyMCMuons>;
using FilteredMyMCMuon = FilteredMyMCMuons::iterator;

using MySmearedElectrons = soa::Join<aod::EMMCParticles, aod::SmearedElectrons>;
using MySmearedElectron = MySmearedElectrons::iterator;

using MySmearedMuons = soa::Join<aod::EMMCParticles, aod::SmearedMuons>;
using MySmearedMuon = MySmearedMuons::iterator;

// template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename... Types>
template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename TLeptons, typename TSmeardMCParitlces>
struct DileptonMC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};
  Configurable<uint> cfgDCAType{"cfgDCAType", 0, "type of DCA for output. 0:3D, 1:XY, 2:Z, else:3D"};
  Configurable<bool> cfgFillUnfolding{"cfgFillUnfolding", false, "flag to fill histograms for unfolding"};
  Configurable<bool> cfgRequireTrueAssociation{"cfgRequireTrueAssociation", false, "flag to require true mc collision association"};
  Configurable<bool> cfgFillSeparateCharmHadronPairs{"cfgFillSeparateCharmHadronPairs", false, "flag to fill different ccbar pairs separately"};

  ConfigurableAxis ConfMllBins{"ConfMllBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00}, "mll bins for output histograms"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTll bins for output histograms"};
  ConfigurableAxis ConfDCAllBins{"ConfDCAllBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAll bins for output histograms"};

  ConfigurableAxis ConfDPtBins{"ConfDPtBins", {220, -1.0, +10.0}, "dpt bins for output histograms"};
  ConfigurableAxis ConfDCAllNarrowBins{"ConfDCAllNarrowBins", {200, 0.0, 10.0}, "narrow DCAll bins for output histograms"};
  ConfigurableAxis ConfTrackDCA{"ConfTrackDCA", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10}, "DCA binning for single tacks"};

  ConfigurableAxis ConfYllBins{"ConfYllBins", {1, -1.f, +1.f}, "yll bins for output histograms"};

  // ConfigurableAxis ConfMmumuBins{"ConfMmumuBins", {VARIABLE_WIDTH, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70, 5.80, 5.90, 6.00, 6.10, 6.20, 6.30, 6.40, 6.50, 6.60, 6.70, 6.80, 6.90, 7.00, 7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.80, 7.90, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.00, 11.50, 12.00}, "mmumu bins for output histograms"}; // for dimuon. one can copy bins here to hyperloop page.

  Configurable<int> cfg_nbin_dphi_ee{"cfg_nbin_dphi_ee", 1, "number of bins for dphi_ee"}; // 36
  Configurable<int> cfg_nbin_deta_ee{"cfg_nbin_deta_ee", 1, "number of bins for deta_ee"}; // 40
  // Configurable<int> cfg_nbin_cos_theta_cs{"cfg_nbin_cos_theta_cs", 1, "number of bins for cos theta cs"}; // 10
  // Configurable<int> cfg_nbin_phi_cs{"cfg_nbin_phi_cs", 1, "number of bins for phi cs"};                   // 10
  Configurable<int> cfg_nbin_aco{"cfg_nbin_aco", 1, "number of bins for acoplanarity"};          // 10
  Configurable<int> cfg_nbin_asym_pt{"cfg_nbin_asym_pt", 1, "number of bins for pt asymmetry"};  // 10
  Configurable<int> cfg_nbin_dphi_e_ee{"cfg_nbin_dphi_e_ee", 1, "number of bins for dphi_ee_e"}; // 18
  ConfigurableAxis ConfPolarizationCosThetaBins{"ConfPolarizationCosThetaBins", {1, -1.f, 1.f}, "cos(theta) bins for polarization analysis"};
  ConfigurableAxis ConfPolarizationPhiBins{"ConfPolarizationPhiBins", {1, -M_PI, M_PI}, "phi bins for polarization analysis"};
  Configurable<int> cfgPolarizationFrame{"cfgPolarizationFrame", 0, "frame of polarization. 0:CS, 1:HX, else:FATAL"};
  ConfigurableAxis ConfPolarizationQuadMomBins{"ConfPolarizationQuadMomBins", {1, -0.5, 1}, "quadrupole moment bins for polarization analysis"}; // quardrupole moment <(3 x cos^2(theta) -1)/2>

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
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

    Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
    Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};
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

    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply phiv cut inherited from prefilter"};
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
    Configurable<bool> rejectWrongMatch{"rejectWrongMatch", false, "flag to reject wrong match between MFT and MCH-MID"}; // this is only for MC study, as we don't know correct match in data.
  } dimuoncuts;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int mRunNumber;
  float d_bz;

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt{"min_mcPt", 0.1, "min. MC pT for generated single lepton"};
    Configurable<float> max_mcPt{"max_mcPt", 1e+10, "max. MC pT for generated single lepton"};
    Configurable<float> min_mcEta{"min_mcEta", -0.8, "max. MC eta for generated single lepton"};
    Configurable<float> max_mcEta{"max_mcEta", +0.8, "max. MC eta for generated single lepton"};
  } mctrackcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view pair_sign_types[3] = {"uls/", "lspp/", "lsmm/"};
  static constexpr std::string_view dilepton_source_types[20] = {
    "sm/Photon/",             // 0
    "sm/PromptPi0/",          // 1
    "sm/NonPromptPi0/",       // 2
    "sm/Eta/",                // 3
    "sm/EtaPrime/",           // 4
    "sm/Rho/",                // 5
    "sm/Omega/",              // 6
    "sm/Phi/",                // 7
    "sm/PromptJPsi/",         // 8
    "sm/NonPromptJPsi/",      // 9
    "sm/PromptPsi2S/",        // 10
    "sm/NonPromptPsi2S/",     // 11
    "sm/Upsilon1S/",          // 12
    "sm/Upsilon2S/",          // 13
    "sm/Upsilon3S/",          // 14
    "ccbar/c2l_c2l/",         // 15
    "bbbar/b2l_b2l/",         // 16
    "bbbar/b2c2l_b2c2l/",     // 17
    "bbbar/b2c2l_b2l_sameb/", // 18
    "bbbar/b2c2l_b2l_diffb/"  // 19
  }; // unordered_map is better, but cannot be constexpr.
  static constexpr std::string_view unfolding_dilepton_source_types[3] = {"sm/", "ccbar/", "bbbar/"};

  ~DileptonMC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms(&fRegistry);
    fRegistry.add("MCEvent/before/hZvtx", "mc vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("MCEvent/before/hZvtx_rec", "rec. mc vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.addClone("MCEvent/before/", "MCEvent/after/");

    std::string mass_axis_title = "m_{ll} (GeV/c^{2})";
    std::string pair_pt_axis_title = "p_{T,ll} (GeV/c)";
    std::string pair_y_axis_title = "y_{ll}";
    std::string pair_dca_axis_title = "DCA_{ll} (#sigma)";
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
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      mass_axis_title = "m_{#mu#mu} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,#mu#mu} (GeV/c)";
      pair_y_axis_title = "y_{#mu#mu}";
      pair_dca_axis_title = "DCA_{#mu#mu}^{XY} (#sigma)";
    }

    // pair info
    const AxisSpec axis_mass{ConfMllBins, mass_axis_title};
    const AxisSpec axis_pt{ConfPtllBins, pair_pt_axis_title};
    const AxisSpec axis_y{ConfYllBins, pair_y_axis_title};
    const AxisSpec axis_dca{ConfDCAllBins, pair_dca_axis_title};
    const AxisSpec axis_pt_meson{ConfPtllBins, "p_{T}^{VM} (GeV/c)"}; // for omega, phi meson pT spectra
    const AxisSpec axis_y_meson{ConfYllBins, "y^{VM}"};               // for omega, phi meson pT spectra

    const AxisSpec axis_dca_narrow{ConfDCAllNarrowBins, pair_dca_axis_title};
    const AxisSpec axis_dpt{ConfDPtBins, "#Delta p_{T,1}^{gen-rec} + #Delta p_{T,2}^{gen-rec} (GeV/c)"};
    const AxisSpec axis_dca_track1{ConfTrackDCA, "DCA_{e,1}^{Z} (#sigma)"};
    const AxisSpec axis_dca_track2{ConfTrackDCA, "DCA_{e,2}^{Z} (#sigma)"};

    std::string frameName = "CS";
    if (cfgPolarizationFrame == 0) {
      frameName = "CS";
    } else if (cfgPolarizationFrame == 1) {
      frameName = "HX";
    } else {
      LOG(fatal) << "set 0 or 1 to cfgPolarizationFrame!";
    }

    const AxisSpec axis_dphi_ee{cfg_nbin_dphi_ee, -M_PI / 2., 3. / 2. * M_PI, "#Delta#varphi = #varphi_{l1} - #varphi_{l2} (rad.)"}; // for kHFll
    const AxisSpec axis_deta_ee{cfg_nbin_deta_ee, -2., 2., "#Delta#eta = #eta_{l1} - #eta_{l2}"};                                    // for kHFll
    const AxisSpec axis_cos_theta_pol{ConfPolarizationCosThetaBins, Form("cos(#theta^{%s})", frameName.data())};                     // for kPolarization, kUPC
    const AxisSpec axis_phi_pol{ConfPolarizationPhiBins, Form("#varphi^{%s} (rad.)", frameName.data())};                             // for kPolarization
    const AxisSpec axis_quadmom{ConfPolarizationQuadMomBins, Form("#frac{3 cos^{2}(#theta^{%s}) -1}{2}", frameName.data())};         // for kPolarization
    const AxisSpec axis_aco{cfg_nbin_aco, 0, 1.f, "#alpha = 1 - #frac{|#varphi_{l^{+}} - #varphi_{l^{-}}|}{#pi}"};                   // for kUPC
    const AxisSpec axis_asym_pt{cfg_nbin_asym_pt, 0, 1.f, "A = #frac{|p_{T,l^{+}} - p_{T,l^{-}}|}{|p_{T,l^{+}} + p_{T,l^{-}}|}"};    // for kUPC
    const AxisSpec axis_dphi_e_ee{cfg_nbin_dphi_e_ee, 0, M_PI, "#Delta#varphi = #varphi_{l} - #varphi_{ll} (rad.)"};                 // for kUPC

    // generated info
    fRegistry.add("Generated/sm/PromptPi0/uls/hs", "gen. dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_y, axis_dphi_ee, axis_deta_ee, axis_cos_theta_pol, axis_phi_pol, axis_quadmom, axis_aco, axis_asym_pt, axis_dphi_e_ee}, true);
    fRegistry.addClone("Generated/sm/PromptPi0/uls/", "Generated/sm/PromptPi0/lspp/");
    fRegistry.addClone("Generated/sm/PromptPi0/uls/", "Generated/sm/PromptPi0/lsmm/");

    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/NonPromptPi0/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Eta/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/EtaPrime/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Rho/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Omega/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Omega2ll/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Phi/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Phi2ll/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/PromptJPsi/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/NonPromptJPsi/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/PromptPsi2S/");
    fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/NonPromptPsi2S/");
    // fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Upsilon1S/");
    // fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Upsilon2S/");
    // fRegistry.addClone("Generated/sm/PromptPi0/", "Generated/sm/Upsilon3S/");

    fRegistry.add("Generated/ccbar/c2l_c2l/uls/hs", "generated dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_y, axis_dphi_ee, axis_deta_ee, axis_cos_theta_pol, axis_phi_pol, axis_quadmom, axis_aco, axis_asym_pt, axis_dphi_e_ee}, true);
    fRegistry.addClone("Generated/ccbar/c2l_c2l/uls/", "Generated/ccbar/c2l_c2l/lspp/");
    fRegistry.addClone("Generated/ccbar/c2l_c2l/uls/", "Generated/ccbar/c2l_c2l/lsmm/");

    fRegistry.addClone("Generated/ccbar/c2l_c2l/", "Generated/bbbar/b2l_b2l/");
    fRegistry.addClone("Generated/ccbar/c2l_c2l/", "Generated/bbbar/b2c2l_b2c2l/");
    fRegistry.addClone("Generated/ccbar/c2l_c2l/", "Generated/bbbar/b2c2l_b2l_sameb/");
    fRegistry.addClone("Generated/ccbar/c2l_c2l/", "Generated/bbbar/b2c2l_b2l_diffb/"); // LS

    // for charmed hadrons // create 28 combinations
    static constexpr std::string_view charmed_mesons[] = {"Dplus", "D0", "Dsplus"}; // 411, 421, 431
    static constexpr std::string_view anti_charmed_mesons[] = {"Dminus", "D0bar", "Dsminus"};
    const int nm = sizeof(charmed_mesons) / sizeof(charmed_mesons[0]);
    static constexpr std::string_view charmed_baryons[] = {"Lcplus", "Xicplus", "Xic0", "Omegac0"}; // 4122, 4232, 4132, 4332
    static constexpr std::string_view anti_charmed_baryons[] = {"Lcminus", "Xicminus", "Xic0bar", "Omegac0bar"};
    const int nb = sizeof(charmed_baryons) / sizeof(charmed_baryons[0]);
    static constexpr std::string_view sum_charmed_mesons[] = {"Dpm", "D0", "Dspm"};
    static constexpr std::string_view sum_charmed_baryons[] = {"Lcpm", "Xicpm", "Xic0", "Omegac0"};

    if (cfgFillSeparateCharmHadronPairs) {
      for (int im = 0; im < nm; im++) {
        fRegistry.addClone("Generated/ccbar/c2l_c2l/", Form("Generated/ccbar/%s_%s/", charmed_mesons[im].data(), anti_charmed_mesons[im].data()));
      }
      for (int ib = 0; ib < nb; ib++) {
        fRegistry.addClone("Generated/ccbar/c2l_c2l/", Form("Generated/ccbar/%s_%s/", charmed_baryons[ib].data(), anti_charmed_baryons[ib].data()));
      }
      for (int im1 = 0; im1 < nm - 1; im1++) {
        for (int im2 = im1 + 1; im2 < nm; im2++) {
          fRegistry.addClone("Generated/ccbar/c2l_c2l/", Form("Generated/ccbar/%s_%s/", sum_charmed_mesons[im1].data(), sum_charmed_mesons[im2].data()));
        }
      }
      for (int ib1 = 0; ib1 < nb - 1; ib1++) {
        for (int ib2 = ib1 + 1; ib2 < nb; ib2++) {
          fRegistry.addClone("Generated/ccbar/c2l_c2l/", Form("Generated/ccbar/%s_%s/", sum_charmed_baryons[ib1].data(), sum_charmed_baryons[ib2].data()));
        }
      }
      for (int im = 0; im < nm; im++) {
        for (int ib = 0; ib < nb; ib++) {
          fRegistry.addClone("Generated/ccbar/c2l_c2l/", Form("Generated/ccbar/%s_%s/", sum_charmed_mesons[im].data(), sum_charmed_baryons[ib].data()));
        }
      }
    }

    // evaluate acceptance for polarization
    fRegistry.add("Generated/VM/All/Phi/hs", "gen. VM #rightarrow ll", kTHnSparseD, {axis_mass, axis_pt, axis_y, axis_cos_theta_pol, axis_phi_pol, axis_quadmom}, true);
    fRegistry.addClone("Generated/VM/All/Phi/", "Generated/VM/All/Rho/");
    fRegistry.addClone("Generated/VM/All/Phi/", "Generated/VM/All/Omega/");
    fRegistry.addClone("Generated/VM/All/Phi/", "Generated/VM/All/PromptJPsi/");
    fRegistry.addClone("Generated/VM/All/Phi/", "Generated/VM/All/NonPromptJPsi/");
    fRegistry.addClone("Generated/VM/All/", "Generated/VM/Acc/");

    // reconstructed pair info
    fRegistry.add("Pair/sm/Photon/uls/hs", "rec. dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_y, axis_dphi_ee, axis_deta_ee, axis_cos_theta_pol, axis_phi_pol, axis_quadmom, axis_aco, axis_asym_pt, axis_dphi_e_ee, axis_dca}, true);

    fRegistry.addClone("Pair/sm/Photon/uls/", "Pair/sm/Photon/lspp/");
    fRegistry.addClone("Pair/sm/Photon/uls/", "Pair/sm/Photon/lsmm/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/PromptPi0/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/NonPromptPi0/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Eta/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/EtaPrime/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Rho/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Omega/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Omega2ll/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Phi/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Phi2ll/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/PromptJPsi/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/NonPromptJPsi/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/PromptPsi2S/");
    fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/NonPromptPsi2S/");
    // fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Upsilon1S/");
    // fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Upsilon2S/");
    // fRegistry.addClone("Pair/sm/Photon/", "Pair/sm/Upsilon3S/");

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      fRegistry.add("Pair/sm/Photon/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0.0f, 1.0f}}, true);
      fRegistry.add("Pair/sm/Photon/uls/hMvsRxy", "m_{ee} vs. r_{xy};r_{xy}^{true} (cm);m_{ee} (GeV/c^{2})", kTH2F, {{100, 0, 100}, {100, 0.0f, 1.0f}}, true);
      for (const auto& strSign : pair_sign_types) {
        fRegistry.add(std::format("Pair/sm/PromptPi0/{0}hMvsPhiV", strSign), "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0.0f, 1.0f}}, true);
        fRegistry.add(std::format("Pair/sm/PromptPi0/{0}hDeltaPtvsDCA", strSign), "#Delta p_{T,1}^{gen-rec} + #Delta p_{T,2}^{gen-rec} vs. DCA_{ee}", kTH2F, {axis_dca_narrow, axis_dpt}, true);
        fRegistry.add(std::format("Pair/sm/PromptPi0/{0}hDCAz1vsDCAz2", strSign), "DCA_{z,1} vs DCA_{z,2}", kTH2F, {axis_dca_track1, axis_dca_track2}, true);

        fRegistry.add(std::format("Pair/sm/NonPromptPi0/{0}hMvsPhiV", strSign), "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0.0f, 1.0f}}, true);
        fRegistry.add(std::format("Pair/sm/NonPromptPi0/{0}hDeltaPtvsDCA", strSign), "#Delta p_{T,1}^{gen-rec} + #Delta p_{T,2}^{gen-rec} vs. DCA_{ee}", kTH2F, {axis_dca_narrow, axis_dpt}, true);
        fRegistry.add(std::format("Pair/sm/NonPromptPi0/{0}hDCAz1vsDCAz2", strSign), "DCA_{z,1} vs DCA_{z,2}", kTH2F, {axis_dca_track1, axis_dca_track2}, true);

        fRegistry.add(std::format("Pair/sm/PromptJPsi/{0}hDeltaPtvsDCA", strSign), "#Delta p_{T,1}^{gen-rec} + #Delta p_{T,2}^{gen-rec} vs. DCA_{ee}", kTH2F, {axis_dca_narrow, axis_dpt}, true);
        fRegistry.add(std::format("Pair/sm/PromptJPsi/{0}hDCAz1vsDCAz2", strSign), "DCA_{z,1} vs DCA_{z,2}", kTH2F, {axis_dca_track1, axis_dca_track2}, true);
        fRegistry.add(std::format("Pair/sm/NonPromptJPsi/{0}hDeltaPtvsDCA", strSign), "#Delta p_{T,1}^{gen-rec} + #Delta p_{T,2}^{gen-rec} vs. DCA_{ee}", kTH2F, {axis_dca_narrow, axis_dpt}, true);
        fRegistry.add(std::format("Pair/sm/NonPromptJPsi/{0}hDCAz1vsDCAz2", strSign), "DCA_{z,1} vs DCA_{z,2}", kTH2F, {axis_dca_track1, axis_dca_track2}, true);
      }
    }

    fRegistry.add("Pair/ccbar/c2l_c2l/uls/hs", "rec. dilepton", kTHnSparseD, {axis_mass, axis_pt, axis_y, axis_dphi_ee, axis_deta_ee, axis_cos_theta_pol, axis_phi_pol, axis_quadmom, axis_aco, axis_asym_pt, axis_dphi_e_ee, axis_dca}, true);
    fRegistry.addClone("Pair/ccbar/c2l_c2l/uls/", "Pair/ccbar/c2l_c2l/lspp/");
    fRegistry.addClone("Pair/ccbar/c2l_c2l/uls/", "Pair/ccbar/c2l_c2l/lsmm/");

    fRegistry.addClone("Pair/ccbar/c2l_c2l/", "Pair/bbbar/b2l_b2l/");
    fRegistry.addClone("Pair/ccbar/c2l_c2l/", "Pair/bbbar/b2c2l_b2c2l/");
    fRegistry.addClone("Pair/ccbar/c2l_c2l/", "Pair/bbbar/b2c2l_b2l_sameb/");
    fRegistry.addClone("Pair/ccbar/c2l_c2l/", "Pair/bbbar/b2c2l_b2l_diffb/"); // LS

    if (cfgFillSeparateCharmHadronPairs) {
      for (int im = 0; im < nm; im++) {
        fRegistry.addClone("Pair/ccbar/c2l_c2l/", Form("Pair/ccbar/%s_%s/", charmed_mesons[im].data(), anti_charmed_mesons[im].data()));
      }
      for (int ib = 0; ib < nb; ib++) {
        fRegistry.addClone("Pair/ccbar/c2l_c2l/", Form("Pair/ccbar/%s_%s/", charmed_baryons[ib].data(), anti_charmed_baryons[ib].data()));
      }
      for (int im1 = 0; im1 < nm - 1; im1++) {
        for (int im2 = im1 + 1; im2 < nm; im2++) {
          fRegistry.addClone("Pair/ccbar/c2l_c2l/", Form("Pair/ccbar/%s_%s/", sum_charmed_mesons[im1].data(), sum_charmed_mesons[im2].data()));
        }
      }
      for (int ib1 = 0; ib1 < nb - 1; ib1++) {
        for (int ib2 = ib1 + 1; ib2 < nb; ib2++) {
          fRegistry.addClone("Pair/ccbar/c2l_c2l/", Form("Pair/ccbar/%s_%s/", sum_charmed_baryons[ib1].data(), sum_charmed_baryons[ib2].data()));
        }
      }
      for (int im = 0; im < nm; im++) {
        for (int ib = 0; ib < nb; ib++) {
          fRegistry.addClone("Pair/ccbar/c2l_c2l/", Form("Pair/ccbar/%s_%s/", sum_charmed_mesons[im].data(), sum_charmed_baryons[ib].data()));
        }
      }
    }

    // for correlated bkg due to mis-identified hadrons, and true combinatorial bkg
    fRegistry.add("Pair/corr_bkg_lh/uls/hs", "rec. bkg", kTHnSparseD, {axis_mass, axis_pt, axis_y, axis_dphi_ee, axis_deta_ee, axis_cos_theta_pol, axis_phi_pol, axis_quadmom, axis_aco, axis_asym_pt, axis_dphi_e_ee, axis_dca}, true);
    fRegistry.addClone("Pair/corr_bkg_lh/uls/", "Pair/corr_bkg_lh/lspp/");
    fRegistry.addClone("Pair/corr_bkg_lh/uls/", "Pair/corr_bkg_lh/lsmm/");
    fRegistry.addClone("Pair/corr_bkg_lh/", "Pair/corr_bkg_hh/");
    fRegistry.addClone("Pair/corr_bkg_lh/", "Pair/comb_bkg/");

    if (doprocessGen_VM) {
      fRegistry.add("Generated/VM/Omega/hPtY", "pT vs. y of #omega", kTH2D, {axis_y_meson, axis_pt_meson}, true); // for pT spectrum
      fRegistry.add("Generated/VM/Phi/hPtY", "pT vs. y of #phi", kTH2D, {axis_y_meson, axis_pt_meson}, true);     // for pT spectrum
    }

    if (cfgFillUnfolding) {
      // for 2D unfolding
      const AxisSpec axis_mass_gen{ConfMllBins, "m_{ll}^{gen} (GeV/c^{2})"};
      const AxisSpec axis_pt_gen{ConfPtllBins, "p_{T,ll}^{gen} (GeV/c)"};
      const AxisSpec axis_mass_rec{ConfMllBins, "m_{ll}^{rec} (GeV/c^{2})"};
      const AxisSpec axis_pt_rec{ConfPtllBins, "p_{T,ll}^{rec} (GeV/c)"};
      fRegistry.add("Unfold/sm/uls/hsRM", "response matrix", kTHnSparseD, {axis_mass_gen, axis_pt_gen, axis_mass_rec, axis_pt_rec}, true);
      fRegistry.add("Unfold/sm/uls/hMiss", "missing dilepton", kTH2D, {axis_mass_gen, axis_pt_gen}, true); // e.g. true eta is in acceptance, but reconstructed eta is out of acceptance.
      fRegistry.add("Unfold/sm/uls/hFake", "fake dilepton", kTH2D, {axis_mass_rec, axis_pt_rec}, true);    // e.g. true eta is out of acceptance, but reconstructed eta is in acceptance.
      fRegistry.addClone("Unfold/sm/uls/", "Unfold/sm/lspp/");
      fRegistry.addClone("Unfold/sm/uls/", "Unfold/sm/lsmm/");
      fRegistry.addClone("Unfold/sm/", "Unfold/ccbar/");
      fRegistry.addClone("Unfold/sm/", "Unfold/bbbar/");
    }
  }

  float beamM1 = o2::constants::physics::MassProton; // mass of beam
  float beamM2 = o2::constants::physics::MassProton; // mass of beam
  float beamE1 = 0.f;                                // beam energy
  float beamE2 = 0.f;                                // beam energy
  float beamP1 = 0.f;                                // beam momentum
  float beamP2 = 0.f;                                // beam momentum

  float leptonM1 = 0.f;
  float leptonM2 = 0.f;
  int pdg_lepton = 0;
  void init(InitContext&)
  {
    if (doprocessAnalysis && doprocessAnalysis_Smeared) {
      LOGF(fatal, "Cannot enable doprocessAnalysis and doprocessAnalysis_Smeared at the same time. Please choose one.");
    }

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);

    DefineEMEventCut();
    addhistograms();

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      DefineDielectronCut();
      leptonM1 = o2::constants::physics::MassElectron;
      leptonM2 = o2::constants::physics::MassElectron;
      pdg_lepton = 11;
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      DefineDimuonCut();
      leptonM1 = o2::constants::physics::MassMuon;
      leptonM2 = o2::constants::physics::MassMuon;
      pdg_lepton = 13;
    }
    if (doprocessNorm) {
      fRegistry.addClone("Event/before/hCollisionCounter", "Event/norm/hCollisionCounter");
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
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, +dielectroncuts.cfg_max_eta_track);
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

    // for pair
    fDimuonCut.SetMassRange(dimuoncuts.cfg_min_mass, dimuoncuts.cfg_max_mass);
    fDimuonCut.SetPairPtRange(dimuoncuts.cfg_min_pair_pt, dimuoncuts.cfg_max_pair_pt);
    fDimuonCut.SetPairYRange(dimuoncuts.cfg_min_pair_y, dimuoncuts.cfg_max_pair_y);
    fDimuonCut.SetPairDCAxyRange(dimuoncuts.cfg_min_pair_dcaxy, dimuoncuts.cfg_max_pair_dcaxy); // DCAxy in cm
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

  template <typename TTrack, typename TMCParticles>
  int FindSMULS(TTrack const& t1mc, TTrack const& t2mc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 22, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 111, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 221, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 331, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 113, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 223, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 333, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 443, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 100443, mcparticles)
      // FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 553, mcparticles),
      // FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 100553, mcparticles),
      // FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, pdg_lepton, 200553, mcparticles)
    };
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <typename TTrack, typename TMCParticles>
  int FindSMLSPP(TTrack const& t1mc, TTrack const& t2mc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, -pdg_lepton, 111, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, -pdg_lepton, 221, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, -pdg_lepton, 331, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, -pdg_lepton, 113, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, -pdg_lepton, 223, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -pdg_lepton, -pdg_lepton, 333, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <typename TTrack, typename TMCParticles>
  int FindSMLSMM(TTrack const& t1mc, TTrack const& t2mc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(t1mc, t2mc, pdg_lepton, pdg_lepton, 111, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, pdg_lepton, pdg_lepton, 221, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, pdg_lepton, pdg_lepton, 331, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, pdg_lepton, pdg_lepton, 113, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, pdg_lepton, pdg_lepton, 223, mcparticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, pdg_lepton, pdg_lepton, 333, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <bool isSmeared, typename T>
  bool isInAcceptance(T const& lepton)
  {
    float pt = 0.f, eta = 0.f;
    if constexpr (isSmeared) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        pt = lepton.ptSmeared();
        eta = lepton.etaSmeared();
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          pt = lepton.ptSmeared_sa_muon();
          eta = lepton.etaSmeared_sa_muon();
        } else if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          pt = lepton.ptSmeared_gl_muon();
          eta = lepton.etaSmeared_gl_muon();
        } else {
          pt = lepton.pt();
          eta = lepton.eta();
        }
      }
    } else {
      pt = lepton.pt();
      eta = lepton.eta();
    }

    return isInAcceptance(pt, eta);

    // if ((mctrackcuts.min_mcPt < pt && pt < mctrackcuts.max_mcPt) && (mctrackcuts.min_mcEta < eta && eta < mctrackcuts.max_mcEta)) {
    //   return true;
    // } else {
    //   return false;
    // }
  }

  bool isInAcceptance(const float pt, const float eta)
  {
    if ((mctrackcuts.min_mcPt < pt && pt < mctrackcuts.max_mcPt) && (mctrackcuts.min_mcEta < eta && eta < mctrackcuts.max_mcEta)) {
      return true;
    } else {
      return false;
    }
  }

  template <int sourceId>
  void fillGenHistograms(const int sign1, const int sign2, const int pdgMotherC1, const int pdgMotherC2, const float mass, const float pt, const float rapidity, const float dphi, const float deta, const float cos_thetaPol, const float phiPol, const float quadmom, const float aco, const float asym, const float dphi_e_ee, const float weight)
  {
    if (sign1 * sign2 < 0) { // ULS
      fRegistry.fill(HIST("Generated/") + HIST(dilepton_source_types[sourceId]) + HIST("uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
    } else if (sign1 > 0 && sign2 > 0) { // LS++
      fRegistry.fill(HIST("Generated/") + HIST(dilepton_source_types[sourceId]) + HIST("lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
    } else if (sign1 < 0 && sign2 < 0) { // LS--
      fRegistry.fill(HIST("Generated/") + HIST(dilepton_source_types[sourceId]) + HIST("lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
    }

    if (dilepton_source_types[sourceId].find("ccbar") != std::string_view::npos && cfgFillSeparateCharmHadronPairs) {
      if (std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 411) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dplus_Dminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dplus_Dminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dplus_Dminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if (std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 421) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/D0_D0bar/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/D0_D0bar/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/D0_D0bar/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if (std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 431) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dsplus_Dsminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dsplus_Dsminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dsplus_Dsminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 421) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 421)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dpm_D0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dpm_D0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dpm_D0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 431) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 431)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Dspm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Dspm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Dspm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 431) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 431)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/D0_Dspm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/D0_Dspm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/D0_Dspm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4122) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Lcplus_Lcminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Lcplus_Lcminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Lcplus_Lcminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4232 && std::abs(pdgMotherC2) == 4232) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Xicplus_Xicminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Xicplus_Xicminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Xicplus_Xicminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4132 && std::abs(pdgMotherC2) == 4132) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Xic0_Xic0bar/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Xic0_Xic0bar/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Xic0_Xic0bar/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4332 && std::abs(pdgMotherC2) == 4332) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Omegac0_Omegac0bar/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Omegac0_Omegac0bar/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Omegac0_Omegac0bar/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 4122 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 4122 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 4122 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Lcpm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4232 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 4232 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Xicpm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Xicpm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Xicpm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4232 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 4232 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Xicpm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Xicpm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Xicpm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4132 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 4132 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Xic0_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Xic0_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Xic0_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4122) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4122)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Lcpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Lcpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Lcpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dpm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4122) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4122)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/D0_Lcpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/D0_Lcpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/D0_Lcpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/D0_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/D0_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/D0_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/D0_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/D0_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/D0_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/D0_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/D0_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/D0_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4122) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4122)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Lcpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Lcpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Lcpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Generated/ccbar/Dspm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, weight);
        }
      }
    }
  }

  template <int sourceId>
  void fillRecHistograms(const int sign1, const int sign2, const int pdgMotherC1, const int pdgMotherC2, const float mass, const float pt, const float rapidity, const float dphi, const float deta, const float cos_thetaPol, const float phiPol, const float quadmom, const float aco, const float asym, const float dphi_e_ee, const float pair_dca, const float weight)
  {
    if (sign1 * sign2 < 0) { // ULS
      fRegistry.fill(HIST("Pair/") + HIST(dilepton_source_types[sourceId]) + HIST("uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
    } else if (sign1 > 0 && sign2 > 0) { // LS++
      fRegistry.fill(HIST("Pair/") + HIST(dilepton_source_types[sourceId]) + HIST("lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
    } else if (sign1 < 0 && sign2 < 0) { // LS--
      fRegistry.fill(HIST("Pair/") + HIST(dilepton_source_types[sourceId]) + HIST("lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
    }

    if (dilepton_source_types[sourceId].find("ccbar") != std::string_view::npos && cfgFillSeparateCharmHadronPairs) {
      if (std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 411) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dplus_Dminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dplus_Dminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dplus_Dminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if (std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 421) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/D0_D0bar/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/D0_D0bar/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/D0_D0bar/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if (std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 431) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dsplus_Dsminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dsplus_Dsminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dsplus_Dsminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 421) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 421)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dpm_D0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dpm_D0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dpm_D0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 431) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 431)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Dspm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Dspm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Dspm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 431) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 431)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/D0_Dspm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/D0_Dspm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/D0_Dspm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4122) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Lcplus_Lcminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Lcplus_Lcminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Lcplus_Lcminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4232 && std::abs(pdgMotherC2) == 4232) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Xicplus_Xicminus/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Xicplus_Xicminus/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Xicplus_Xicminus/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4132 && std::abs(pdgMotherC2) == 4132) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Xic0_Xic0bar/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Xic0_Xic0bar/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Xic0_Xic0bar/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if (std::abs(pdgMotherC1) == 4332 && std::abs(pdgMotherC2) == 4332) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Omegac0_Omegac0bar/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Omegac0_Omegac0bar/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Omegac0_Omegac0bar/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 4122 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 4122 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4122 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 4122 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Lcpm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4232 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 4232 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Xicpm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Xicpm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Xicpm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4232 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 4232 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Xicpm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Xicpm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Xicpm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 4132 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 4132 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Xic0_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Xic0_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Xic0_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4122) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4122)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Lcpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Lcpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Lcpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 411 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 411 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dpm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4122) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4122)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/D0_Lcpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/D0_Lcpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/D0_Lcpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/D0_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/D0_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/D0_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/D0_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/D0_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/D0_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 421 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 421 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/D0_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/D0_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/D0_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4122) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4122)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Lcpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Lcpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Lcpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4232) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4232)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Xicpm/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Xicpm/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Xicpm/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4132) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4132)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Xic0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Xic0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Xic0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      } else if ((std::abs(pdgMotherC1) == 431 && std::abs(pdgMotherC2) == 4332) || (std::abs(pdgMotherC2) == 431 && std::abs(pdgMotherC1) == 4332)) {
        if (sign1 * sign2 < 0) { // ULS
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Omegac0/uls/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 > 0 && sign2 > 0) { // LS++
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Omegac0/lspp/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        } else if (sign1 < 0 && sign2 < 0) { // LS--
          fRegistry.fill(HIST("Pair/ccbar/Dspm_Omegac0/lsmm/hs"), mass, pt, rapidity, dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, dphi_e_ee, pair_dca, weight);
        }
      }
    }
  }

  template <bool isSmeared, typename TCollision, typename TMCCollisions, typename TTrack1, typename TTrack2, typename TCut, typename TAllTracks, typename TMCParticles>
  bool fillTruePairInfo(TCollision const& collision, TMCCollisions const&, TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TAllTracks const& tracks, TMCParticles const& mcparticles)
  {
    auto t1mc = mcparticles.iteratorAt(t1.emmcparticleId());
    auto t2mc = mcparticles.iteratorAt(t2.emmcparticleId());
    bool is_pair_from_same_mcevent = (t1mc.emmceventId() == t2mc.emmceventId());

    auto mccollision1 = t1mc.template emmcevent_as<TMCCollisions>();
    auto mccollision2 = t2mc.template emmcevent_as<TMCCollisions>();
    if (cfgEventGeneratorType >= 0 && mccollision1.getSubGeneratorId() != cfgEventGeneratorType) {
      return false;
    }
    if (cfgEventGeneratorType >= 0 && mccollision2.getSubGeneratorId() != cfgEventGeneratorType) {
      return false;
    }
    if (!isInAcceptance<isSmeared>(t1mc) || !isInAcceptance<isSmeared>(t2mc)) {
      return false;
    }

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
      if (!cut.IsSelectedPair(t1, t2, d_bz, 0.0)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedTrack<false>(t1) || !cut.template IsSelectedTrack<false>(t2)) {
        return false;
      }

      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch<false>(t1, cut, tracks)) {
        return false;
      }
      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch<false>(t2, cut, tracks)) {
        return false;
      }
      if (dimuoncuts.rejectWrongMatch) {
        if (t1.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && t1.emmcparticleId() != t1.emmftmcparticleId()) {
          return false;
        }
        if (t2.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && t2.emmcparticleId() != t2.emmftmcparticleId()) {
          return false;
        }
      }

      if (!cut.IsSelectedPair(t1, t2)) {
        return false;
      }
    }

    float pt1 = 0.f, eta1 = 0.f, phi1 = 0.f, pt2 = 0.f, eta2 = 0.f, phi2 = 0.f;
    if constexpr (isSmeared) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        pt1 = t1mc.ptSmeared();
        eta1 = t1mc.etaSmeared();
        phi1 = t1mc.phiSmeared();
        pt2 = t2mc.ptSmeared();
        eta2 = t2mc.etaSmeared();
        phi2 = t2mc.phiSmeared();
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          pt1 = t1mc.ptSmeared_sa_muon();
          eta1 = t1mc.etaSmeared_sa_muon();
          phi1 = t1mc.phiSmeared_sa_muon();
          pt2 = t2mc.ptSmeared_sa_muon();
          eta2 = t2mc.etaSmeared_sa_muon();
          phi2 = t2mc.phiSmeared_sa_muon();
        } else if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          pt1 = t1mc.ptSmeared_gl_muon();
          eta1 = t1mc.etaSmeared_gl_muon();
          phi1 = t1mc.phiSmeared_gl_muon();
          pt2 = t2mc.ptSmeared_gl_muon();
          eta2 = t2mc.etaSmeared_gl_muon();
          phi2 = t2mc.phiSmeared_gl_muon();
        } else {
          pt1 = t1mc.pt();
          eta1 = t1mc.eta();
          phi1 = t1mc.phi();
          pt2 = t2mc.pt();
          eta2 = t2mc.eta();
          phi2 = t2mc.phi();
        }
      }
    } else {
      pt1 = t1mc.pt();
      eta1 = t1mc.eta();
      phi1 = t1mc.phi();
      pt2 = t2mc.pt();
      eta2 = t2mc.eta();
      phi2 = t2mc.phi();
    }

    ROOT::Math::PtEtaPhiMVector v1mc(pt1, eta1, phi1, leptonM1);
    ROOT::Math::PtEtaPhiMVector v2mc(pt2, eta2, phi2, leptonM2);
    ROOT::Math::PtEtaPhiMVector v12mc = v1mc + v2mc;

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (v12mc.Rapidity() < dielectroncuts.cfg_min_pair_y || dielectroncuts.cfg_max_pair_y < v12mc.Rapidity()) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (v12mc.Rapidity() < dimuoncuts.cfg_min_pair_y || dimuoncuts.cfg_max_pair_y < v12mc.Rapidity()) {
        return false;
      }
    }

    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[std::make_pair(t1.globalIndex(), t2.globalIndex())];
    }
    // LOGF(info, "t1.sign() = %d, t2.sign() = %d, map_weight[std::make_pair(%d, %d)] = %f", t1.sign(), t2.sign(), t1.globalIndex(), t2.globalIndex(), weight);

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
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);

    float deta = v1.Eta() - v2.Eta();
    float dphi = v1.Phi() - v2.Phi();
    o2::math_utils::bringToPMPi(dphi);

    float aco = 1.f - std::fabs(dphi) / M_PI;
    float asym = std::fabs(v1.Pt() - v2.Pt()) / (v1.Pt() + v2.Pt());
    float dphi_e_ee = v1.Phi() - v12.Phi();
    o2::math_utils::bringToPMPi(dphi_e_ee);
    dphi = RecoDecay::constrainAngle(dphi, -o2::constants::math::PIHalf, 1); // shift dphi in [-pi/2, +3pi/2] rad.

    std::array<float, 4> arrP1 = {t1.px(), t1.py(), t1.pz(), leptonM1};
    std::array<float, 4> arrP2 = {t2.px(), t2.py(), t2.pz(), leptonM2};
    float cos_thetaPol = 999, phiPol = 999.f;
    if (cfgPolarizationFrame == 0) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, t1.sign(), cos_thetaPol, phiPol);
    } else if (cfgPolarizationFrame == 1) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, t1.sign(), cos_thetaPol, phiPol);
    }
    o2::math_utils::bringToPMPi(phiPol);
    float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;

    if ((FindCommonMotherFrom2ProngsWithoutPDG(t1mc, t2mc) > 0 || IsHF(t1mc, t2mc, mcparticles) > 0) && is_pair_from_same_mcevent) { // for bkg study
      if (std::abs(t1mc.pdgCode()) != pdg_lepton || std::abs(t2mc.pdgCode()) != pdg_lepton) {                                        // hh or lh correlated bkg
        if (std::abs(t1mc.pdgCode()) != pdg_lepton && std::abs(t2mc.pdgCode()) != pdg_lepton) {                                      // hh correlated bkg
          if (t1.sign() * t2.sign() < 0) {                                                                                           // ULS
            fRegistry.fill(HIST("Pair/corr_bkg_hh/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
          } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
            fRegistry.fill(HIST("Pair/corr_bkg_hh/lspp/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
          } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
            fRegistry.fill(HIST("Pair/corr_bkg_hh/lsmm/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
          }
        } else {                           // lh correlated bkg
          if (t1.sign() * t2.sign() < 0) { // ULS
            fRegistry.fill(HIST("Pair/corr_bkg_lh/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
          } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
            fRegistry.fill(HIST("Pair/corr_bkg_lh/lspp/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
          } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
            fRegistry.fill(HIST("Pair/corr_bkg_lh/lsmm/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
          }
        }
      }
    } else {                           // true combinatorial bkg
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/comb_bkg/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/comb_bkg/lspp/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/comb_bkg/lsmm/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight);
      }
    }

    if (std::abs(t1mc.pdgCode()) != pdg_lepton || std::abs(t2mc.pdgCode()) != pdg_lepton) {
      return false;
    }

    if (!is_pair_from_same_mcevent) {
      return false;
    }
    if (cfgRequireTrueAssociation && (t1mc.emmceventId() != collision.emmceventId() || t2mc.emmceventId() != collision.emmceventId())) {
      return false;
    }
    int mother_id = std::max({FindSMULS(t1mc, t2mc, mcparticles), FindSMULS(t2mc, t1mc, mcparticles), FindSMLSPP(t1mc, t2mc, mcparticles), FindSMLSMM(t1mc, t2mc, mcparticles)});
    int hfee_type = IsHF(t1mc, t2mc, mcparticles);
    if (mother_id < 0 && hfee_type < 0) {
      return false;
    }

    if (mother_id > -1) {
      auto mcmother = mcparticles.iteratorAt(mother_id);
      if (mcmother.isPhysicalPrimary() || mcmother.producedByGenerator()) {
        if ((t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && (t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
          float deltaPt1 = t1mc.pt() - t1.pt();
          float deltaPt2 = t2mc.pt() - t2.pt();
          switch (std::abs(mcmother.pdgCode())) {
            case 111:
              if (IsFromCharm(mcmother, mcparticles) < 0 && IsFromBeauty(mcmother, mcparticles) < 0) {                                                                                             // prompt pi0
                fillRecHistograms<1>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // prompt pi0
                if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
                  if (t1.sign() * t2.sign() < 0) { // ULS
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/uls/hMvsPhiV"), phiv, v12.M());
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/uls/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/uls/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/lspp/hMvsPhiV"), phiv, v12.M());
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/lspp/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/lspp/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/lsmm/hMvsPhiV"), phiv, v12.M());
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/lsmm/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/PromptPi0/lsmm/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  }
                }
              } else {                                                                                                                                                                             // non-prompt pi0
                fillRecHistograms<2>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // non-prompt pi0
                if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
                  if (t1.sign() * t2.sign() < 0) { // ULS
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/uls/hMvsPhiV"), phiv, v12.M());
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/uls/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/uls/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/lspp/hMvsPhiV"), phiv, v12.M());
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/lspp/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/lspp/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/lsmm/hMvsPhiV"), phiv, v12.M());
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/lsmm/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/NonPromptPi0/lsmm/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  }
                }
              }
              break;
            case 221:
              fillRecHistograms<3>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // eta
              break;
            case 331:
              fillRecHistograms<4>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // eta'
              break;
            case 113:
              fillRecHistograms<5>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // rho
              break;
            case 223:
              fillRecHistograms<6>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // omega
              if (mcmother.daughtersIds().size() == 2) {
                if (t1.sign() * t2.sign() < 0) {                                                                                                                                                    // ULS
                  fRegistry.fill(HIST("Pair/sm/Omega2ll/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // omeag->ee
                }
              }
              break;
            case 333:
              fillRecHistograms<7>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // phi
              if (mcmother.daughtersIds().size() == 2) {
                if (t1.sign() * t2.sign() < 0) {                                                                                                                                                  // ULS
                  fRegistry.fill(HIST("Pair/sm/Phi2ll/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // phi->ee
                }
              }
              break;
            case 443:
              if (IsFromBeauty(mcmother, mcparticles) > 0) {
                fillRecHistograms<9>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // non-prompt J/psi
                if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
                  if (t1.sign() * t2.sign() < 0) { // ULS
                    fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/uls/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/uls/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
                    fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/lspp/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/lspp/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
                    fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/lsmm/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/lsmm/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  }
                }
              } else {
                fillRecHistograms<8>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // prompt J/psi
                if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
                  if (t1.sign() * t2.sign() < 0) { // ULS
                    fRegistry.fill(HIST("Pair/sm/PromptJPsi/uls/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/PromptJPsi/uls/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
                    fRegistry.fill(HIST("Pair/sm/PromptJPsi/lspp/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/PromptJPsi/lspp/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
                    fRegistry.fill(HIST("Pair/sm/PromptJPsi/lsmm/hDeltaPtvsDCA"), pair_dca, deltaPt1 + deltaPt2);
                    fRegistry.fill(HIST("Pair/sm/PromptJPsi/lsmm/hDCAz1vsDCAz2"), dcaZinSigma(t1), dcaZinSigma(t2));
                  }
                }
              }
              break;
            case 100443:
              if (IsFromBeauty(mcmother, mcparticles) > 0) {
                fillRecHistograms<11>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // non-prompt psi2S
              } else {
                fillRecHistograms<10>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // prompt psi2S
              }
              break;
            default:
              break;
          }
        } else if (!(t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && !(t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
          switch (std::abs(mcmother.pdgCode())) {
            case 22:
              fillRecHistograms<0>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // photon conversion
              if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
                fRegistry.fill(HIST("Pair/sm/Photon/uls/hMvsPhiV"), phiv, v12.M());
                float rxy_gen = std::sqrt(std::pow(t1mc.vx(), 2) + std::pow(t1mc.vy(), 2));
                fRegistry.fill(HIST("Pair/sm/Photon/uls/hMvsRxy"), rxy_gen, v12.M());
              }
              break;
            default:
              break;
          }
        } // end of primary/secondary selection
      } // end of primary selection for same mother
    } else if (hfee_type > -1) {
      if ((t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && (t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
        auto mp1 = mcparticles.iteratorAt(t1mc.mothersIds()[0]);
        auto mp2 = mcparticles.iteratorAt(t2mc.mothersIds()[0]);
        switch (hfee_type) {
          case static_cast<int>(EM_HFeeType::kCe_Ce):
            fillRecHistograms<15>(t1.sign(), t2.sign(), mp1.pdgCode(), mp2.pdgCode(), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // c2l_c2l
            break;
          case static_cast<int>(EM_HFeeType::kBe_Be):
            fillRecHistograms<16>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // b2l_b2l
            break;
          case static_cast<int>(EM_HFeeType::kBCe_BCe):
            fillRecHistograms<17>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // b2c2l_b2c2l
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
            fillRecHistograms<18>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // b2c2l_b2l_sameb
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
            fillRecHistograms<19>(t1.sign(), t2.sign(), 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), pair_dca, weight); // b2c2l_b2l_diffb
            break;
          default:
            break;
        }
      }
    } // end of HF evaluation
    return true;
  }

  template <bool isSmeared, typename TMCParticle, typename TMCParticles>
  bool fillGenPairInfo(TMCParticle const& t1, TMCParticle const& t2, TMCParticles const& mcparticles)
  {
    if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
      return false;
    }
    if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
      return false;
    }

    int mother_id = std::max({FindSMULS(t1, t2, mcparticles), FindSMULS(t2, t1, mcparticles), FindSMLSPP(t1, t2, mcparticles), FindSMLSMM(t1, t2, mcparticles)});
    int hfee_type = IsHF(t1, t2, mcparticles);
    if (mother_id < 0 && hfee_type < 0) {
      return false;
    }

    float pt1 = 0.f, eta1 = 0.f, phi1 = 0.f, pt2 = 0.f, eta2 = 0.f, phi2 = 0.f;
    if constexpr (isSmeared) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        pt1 = t1.ptSmeared();
        eta1 = t1.etaSmeared();
        phi1 = t1.phiSmeared();
        pt2 = t2.ptSmeared();
        eta2 = t2.etaSmeared();
        phi2 = t2.phiSmeared();
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          pt1 = t1.ptSmeared_sa_muon();
          eta1 = t1.etaSmeared_sa_muon();
          phi1 = t1.phiSmeared_sa_muon();
          pt2 = t2.ptSmeared_sa_muon();
          eta2 = t2.etaSmeared_sa_muon();
          phi2 = t2.phiSmeared_sa_muon();
        } else if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          pt1 = t1.ptSmeared_gl_muon();
          eta1 = t1.etaSmeared_gl_muon();
          phi1 = t1.phiSmeared_gl_muon();
          pt2 = t2.ptSmeared_gl_muon();
          eta2 = t2.etaSmeared_gl_muon();
          phi2 = t2.phiSmeared_gl_muon();
        } else {
          pt1 = t1.pt();
          eta1 = t1.eta();
          phi1 = t1.phi();
          pt2 = t2.pt();
          eta2 = t2.eta();
          phi2 = t2.phi();
        }
      }
    } else {
      pt1 = t1.pt();
      eta1 = t1.eta();
      phi1 = t1.phi();
      pt2 = t2.pt();
      eta2 = t2.eta();
      phi2 = t2.phi();
    }

    ROOT::Math::PtEtaPhiMVector v1(pt1, eta1, phi1, leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(pt2, eta2, phi2, leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float deta = v1.Eta() - v2.Eta();
    float dphi = v1.Phi() - v2.Phi();
    o2::math_utils::bringToPMPi(dphi);

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (v12.Rapidity() < dielectroncuts.cfg_min_pair_y || dielectroncuts.cfg_max_pair_y < v12.Rapidity()) {
        return false;
      }
      // if (dielectroncuts.cfg_apply_detadphi && std::pow(deta / dielectroncuts.cfg_min_deta, 2) + std::pow(dphi / dielectroncuts.cfg_min_dphi, 2) < 1.f) {
      //   continue;
      // }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (v12.Rapidity() < dimuoncuts.cfg_min_pair_y || dimuoncuts.cfg_max_pair_y < v12.Rapidity()) {
        return false;
      }
      // if (dimuoncuts.cfg_apply_detadphi && std::pow(deta / dimuoncuts.cfg_min_deta, 2) + std::pow(dphi / dimuoncuts.cfg_min_dphi, 2) < 1.f) {
      //   continue;
      // }
    }

    int sign1 = -t1.pdgCode() / pdg_lepton;
    int sign2 = -t2.pdgCode() / pdg_lepton;

    float aco = 1.f - std::fabs(dphi) / M_PI;
    float asym = std::fabs(v1.Pt() - v2.Pt()) / (v1.Pt() + v2.Pt());
    float dphi_e_ee = v1.Phi() - v12.Phi();
    o2::math_utils::bringToPMPi(dphi_e_ee);
    dphi = RecoDecay::constrainAngle(dphi, -o2::constants::math::PIHalf, 1); // shift dphi in [-pi/2, +3pi/2] rad. after deta-dphi cut.

    std::array<float, 4> arrP1 = {static_cast<float>(v1.Px()), static_cast<float>(v1.Py()), static_cast<float>(v1.Pz()), leptonM1};
    std::array<float, 4> arrP2 = {static_cast<float>(v2.Px()), static_cast<float>(v2.Py()), static_cast<float>(v2.Pz()), leptonM2};
    float cos_thetaPol = 999, phiPol = 999.f;
    if (cfgPolarizationFrame == 0) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, -t1.pdgCode() / pdg_lepton, cos_thetaPol, phiPol);
    } else if (cfgPolarizationFrame == 1) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, -t1.pdgCode() / pdg_lepton, cos_thetaPol, phiPol);
    }
    o2::math_utils::bringToPMPi(phiPol);
    float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;

    if (!isInAcceptance<isSmeared>(t1) || !isInAcceptance<isSmeared>(t2)) {
      return false;
    }

    float weight = 1.f;
    if (mother_id > -1) {
      auto mcmother = mcparticles.iteratorAt(mother_id);
      if (mcmother.isPhysicalPrimary() || mcmother.producedByGenerator()) {
        switch (std::abs(mcmother.pdgCode())) {
          case 111:
            if (IsFromCharm(mcmother, mcparticles) < 0 && IsFromBeauty(mcmother, mcparticles) < 0) {                                                                           // prompt pi0
              fillGenHistograms<1>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // prompt pi0
            } else {                                                                                                                                                           // non-prompt pi0
              fillGenHistograms<2>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // non-prompt pi0
            }
            break;
          case 221:
            fillGenHistograms<3>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // eta
            break;
          case 331:
            fillGenHistograms<4>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // eta'
            break;
          case 113:
            fillGenHistograms<5>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // rho
            break;
          case 223:
            fillGenHistograms<6>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // omega
            if (mcmother.daughtersIds().size() == 2) {
              fRegistry.fill(HIST("Generated/sm/Omega2ll/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // omega->ee
            }
            break;
          case 333:
            fillGenHistograms<7>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // phi
            if (mcmother.daughtersIds().size() == 2) {
              fRegistry.fill(HIST("Generated/sm/Phi2ll/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee)); // phi->ee
            }
            break;
          case 443:
            if (IsFromBeauty(mcmother, mcparticles) > 0) {
              fillGenHistograms<9>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // non-prompt J/psi
            } else {
              fillGenHistograms<8>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // prompt J/psi
            }
            break;
          case 100443:
            if (IsFromBeauty(mcmother, mcparticles) > 0) {
              fillGenHistograms<11>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // non-prompt psi2S
            } else {
              fillGenHistograms<10>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // prompt psi2S
            }
            break;
          default:
            break;
        }
      }
    } else if (hfee_type > -1) {
      auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
      auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
      switch (hfee_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce):
          fillGenHistograms<15>(sign1, sign2, mp1.pdgCode(), mp2.pdgCode(), v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // c2l_c2l
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be):
          fillGenHistograms<16>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // b2l_b2l
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe):
          fillGenHistograms<17>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // b2c2l_b2c2l
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
          fillGenHistograms<18>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // b2c2l_b2l_sameb
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
          fillGenHistograms<19>(sign1, sign2, 0, 0, v12.M(), v12.Pt(), v12.Rapidity(), dphi, deta, cos_thetaPol, phiPol, quadmom, aco, asym, std::fabs(dphi_e_ee), weight); // b2c2l_b2l_diffb
          break;
        default:
          break;
      }
    } // end of HF evaluation
    return false;
  }

  template <typename TMCParticle, typename TMCParticles>
  bool fillGenParticleAcc(TMCParticle const& mcParticle, TMCParticles const& mcParticles)
  {
    if (!mcParticle.isPhysicalPrimary() && !mcParticle.producedByGenerator()) {
      return false;
    }
    if (mcParticle.daughtersIds().size() < 2) {
      return false;
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (mcParticle.y() < dielectroncuts.cfg_min_pair_y || dielectroncuts.cfg_max_pair_y < mcParticle.y()) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (mcParticle.y() < dimuoncuts.cfg_min_pair_y || dimuoncuts.cfg_max_pair_y < mcParticle.y()) {
        return false;
      }
    }

    int pdg = mcParticle.pdgCode();
    if (std::abs(pdg) == 113 && mcParticle.daughtersIds().size() != 2) { // reject dalitz decay
      return false;
    }
    if (std::abs(pdg) == 223 && mcParticle.daughtersIds().size() != 2) { // reject dalitz decay
      return false;
    }
    if (std::abs(pdg) == 333 && mcParticle.daughtersIds().size() != 2) { // reject dalitz decay
      return false;
    }
    // accept radiative decay of charmonia (ee + multiple gamma).

    // float pt1 = 0.f, eta1 = 0.f, phi1 = 0.f, sign1 = 0.f;
    // float pt2 = 0.f, eta2 = 0.f, phi2 = 0.f, sign2 = 0.f;
    std::vector<std::array<float, 4>> vDau;
    vDau.reserve(mcParticle.daughtersIds().size());
    for (const auto& daughterId : mcParticle.daughtersIds()) {
      auto dau = mcParticles.rawIteratorAt(daughterId);
      if (std::abs(dau.pdgCode()) == pdg_lepton) {
        vDau.emplace_back(std::array<float, 4>{dau.pt(), dau.eta(), dau.phi(), dau.pdgCode() > 0 ? -1.f : +1.f});
      }
    }
    if (vDau.size() != 2 || vDau[0][3] * vDau[1][3] > 0.f) { // decay objects must be ULS 2 leptons.
      vDau.clear();
      vDau.shrink_to_fit();
      return false;
    }

    // LOGF(info, "mcParticle.globalIndex() = %d, mcParticle.pdgCode() = %d, mcParticle.daughtersIds().size() = %d", mcParticle.globalIndex(), mcParticle.pdgCode(), mcParticle.daughtersIds().size());

    ROOT::Math::PtEtaPhiMVector v1(vDau[0][0], vDau[0][1], vDau[0][2], leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(vDau[1][0], vDau[1][1], vDau[1][2], leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    std::array<float, 4> arrP1 = {static_cast<float>(v1.Px()), static_cast<float>(v1.Py()), static_cast<float>(v1.Pz()), leptonM1};
    std::array<float, 4> arrP2 = {static_cast<float>(v2.Px()), static_cast<float>(v2.Py()), static_cast<float>(v2.Pz()), leptonM2};
    float cos_thetaPol = 999, phiPol = 999.f;
    if (cfgPolarizationFrame == 0) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, vDau[0][3] > 0 ? 1 : -1, cos_thetaPol, phiPol);
    } else if (cfgPolarizationFrame == 1) {
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, vDau[0][3] > 0 ? 1 : -1, cos_thetaPol, phiPol);
    }
    o2::math_utils::bringToPMPi(phiPol);
    float quadmom = (3.f * std::pow(cos_thetaPol, 2) - 1.f) / 2.f;

    float weight = 1.f;
    switch (std::abs(mcParticle.pdgCode())) {
      case 113:
        fRegistry.fill(HIST("Generated/VM/All/Rho/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
        if (isInAcceptance(v1.Pt(), v1.Eta()) && isInAcceptance(v2.Pt(), v2.Eta())) {
          fRegistry.fill(HIST("Generated/VM/Acc/Rho/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
        }
        break;
      case 223:
        fRegistry.fill(HIST("Generated/VM/All/Omega/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
        if (isInAcceptance(v1.Pt(), v1.Eta()) && isInAcceptance(v2.Pt(), v2.Eta())) {
          fRegistry.fill(HIST("Generated/VM/Acc/Omega/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
        }
        break;
      case 333:
        fRegistry.fill(HIST("Generated/VM/All/Phi/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
        if (isInAcceptance(v1.Pt(), v1.Eta()) && isInAcceptance(v2.Pt(), v2.Eta())) {
          fRegistry.fill(HIST("Generated/VM/Acc/Phi/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
        }
        break;
      case 443:
        if (IsFromBeauty(mcParticle, mcParticles) > 0) {
          fRegistry.fill(HIST("Generated/VM/All/NonPromptJPsi/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
          if (isInAcceptance(v1.Pt(), v1.Eta()) && isInAcceptance(v2.Pt(), v2.Eta())) {
            fRegistry.fill(HIST("Generated/VM/Acc/NonPromptJPsi/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
          }
        } else {
          fRegistry.fill(HIST("Generated/VM/All/PromptJPsi/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
          if (isInAcceptance(v1.Pt(), v1.Eta()) && isInAcceptance(v2.Pt(), v2.Eta())) {
            fRegistry.fill(HIST("Generated/VM/Acc/PromptJPsi/hs"), v12.M(), mcParticle.pt(), mcParticle.y(), cos_thetaPol, phiPol, quadmom, weight);
          }
        }
        break;
      default:
        break;
    }

    vDau.clear();
    vDau.shrink_to_fit();
    return false;
  }

  template <bool isSmeared, typename TCollisions, typename TMCLeptons, typename TPreslice, typename TCut, typename TAllTracks, typename TMCCollisions, typename TMCParticles>
  void runTruePairing(TCollisions const& collisions, TMCLeptons const& posTracks, TMCLeptons const& negTracks, TPreslice const& perCollision, TCut const& cut, TAllTracks const& tracks, TMCCollisions const& mccollisions, TMCParticles const& mcparticles)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
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

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillTruePairInfo<isSmeared>(collision, mccollisions, pos, neg, cut, tracks, mcparticles);
      } // end of ULS pair loop

      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        fillTruePairInfo<isSmeared>(collision, mccollisions, pos1, pos2, cut, tracks, mcparticles);
      } // end of LS++ pair loop

      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        fillTruePairInfo<isSmeared>(collision, mccollisions, neg1, neg2, cut, tracks, mcparticles);
      } // end of LS-- pair loop

    } // end of collision loop
  }

  template <bool isSmeared, typename TCollisions, typename TMCCollisions, typename TMCLeptons, typename TMCParticles>
  void runGenInfo(TCollisions const& collisions, TMCCollisions const& mccollisions, TMCLeptons const& posTracksMC, TMCLeptons const& negTracksMC, TMCParticles const& mcparticles)
  {
    for (const auto& mccollision : mccollisions) {
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }
      fRegistry.fill(HIST("MCEvent/before/hZvtx"), mccollision.posZ());
      if (mccollision.mpemeventId() < 0) {
        continue;
      }
      auto collision = collisions.rawIteratorAt(mccollision.mpemeventId());

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("MCEvent/before/hZvtx_rec"), mccollision.posZ());
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      fRegistry.fill(HIST("MCEvent/after/hZvtx"), mccollision.posZ());

      auto posTracks_per_coll = posTracksMC.sliceByCachedUnsorted(aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC.sliceByCachedUnsorted(aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);

      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillGenPairInfo<isSmeared>(t1, t2, mcparticles);
      } // end of true ULS pair loop

      for (const auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        fillGenPairInfo<isSmeared>(t1, t2, mcparticles);
      } // end of true LS++ pair loop

      for (const auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        fillGenPairInfo<isSmeared>(t1, t2, mcparticles);
      } // end of true LS-- pair loop

      // acceptance for polarization of vector mesons
      auto mcParticles_per_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (const auto& mcParticle : mcParticles_per_coll) {
        if (!mcParticle.isPhysicalPrimary() && !mcParticle.producedByGenerator()) {
          continue;
        }
        int pdg = std::abs(mcParticle.pdgCode());
        if (pdg == 113 || pdg == 223 || pdg == 333 || pdg == 443) { // select only VMs
          fillGenParticleAcc(mcParticle, mcparticles);              // VMs that decay into dilepton are stored in derived data. This is sufficient for polarization.
        }
      } // end of mc particle loop
    } // end of collision loop
  }

  template <bool is_wo_acc = false, typename TTrack1, typename TTrack2, typename TCut, typename TAllTracks>
  bool isPairOK(TTrack1 const& t1, TTrack2 const& t2, TCut const& cut, TAllTracks const& tracks)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
        if (!cut.template IsSelectedTrack<is_wo_acc>(t1) || !cut.template IsSelectedTrack<is_wo_acc>(t2)) {
          return false;
        }
      } else { // cut-based
        if (!cut.template IsSelectedTrack<is_wo_acc>(t1) || !cut.template IsSelectedTrack<is_wo_acc>(t2)) {
          return false;
        }
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedTrack<is_wo_acc>(t1) || !cut.template IsSelectedTrack<is_wo_acc>(t2)) {
        return false;
      }
      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch<is_wo_acc>(t1, cut, tracks)) {
        return false;
      }
      if (!o2::aod::pwgem::dilepton::utils::emtrackutil::isBestMatch<is_wo_acc>(t2, cut, tracks)) {
        return false;
      }
    }

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (!cut.template IsSelectedPair<is_wo_acc>(t1, t2, d_bz, 0.0)) {
        return false;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (!cut.template IsSelectedPair<is_wo_acc>(t1, t2)) {
        return false;
      }
    }
    return true;
  }

  std::map<std::pair<int, int>, float> map_weight; // <posId, negId> -> float
  template <typename TCollisions, typename TTracks1, typename TTracks2, typename TPresilce, typename TCut, typename TAllTracks, typename TMCCollisions, typename TMCParticles>
  void fillPairWeightMap(TCollisions const& collisions, TTracks1 const& posTracks, TTracks2 const& negTracks, TPresilce const& perCollision, TCut const& cut, TAllTracks const& tracks, TMCCollisions const&, TMCParticles const& mcparticles)
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

      // auto mccollision = collision.template emmcevent_as<TMCCollisions>();
      // if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
      //   continue;
      // }

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache);

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        auto mcpos = mcparticles.iteratorAt(pos.emmcparticleId());
        auto mccollision_from_pos = mcpos.template emmcevent_as<TMCCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_pos.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        auto mcneg = mcparticles.iteratorAt(neg.emmcparticleId());
        auto mccollision_from_neg = mcneg.template emmcevent_as<TMCCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_neg.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }

        if (isPairOK(pos, neg, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(pos.globalIndex(), neg.globalIndex()));
        }
      }
      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        auto mcpos1 = mcparticles.iteratorAt(pos1.emmcparticleId());
        auto mccollision_from_pos1 = mcpos1.template emmcevent_as<TMCCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_pos1.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        auto mcpos2 = mcparticles.iteratorAt(pos2.emmcparticleId());
        auto mccollision_from_pos2 = mcpos2.template emmcevent_as<TMCCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_pos2.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }

        if (isPairOK(pos1, pos2, cut, tracks)) {
          passed_pairIds.emplace_back(std::make_pair(pos1.globalIndex(), pos2.globalIndex()));
        }
      }
      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        auto mcneg1 = mcparticles.iteratorAt(neg1.emmcparticleId());
        auto mccollision_from_neg1 = mcneg1.template emmcevent_as<TMCCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_neg1.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
        auto mcneg2 = mcparticles.iteratorAt(neg2.emmcparticleId());
        auto mccollision_from_neg2 = mcneg2.template emmcevent_as<TMCCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_neg2.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }
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

  template <typename TTrack>
  bool isPairInAcc(TTrack const& t1, TTrack const& t2)
  {
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if ((t1.pt() < dielectroncuts.cfg_min_pt_track || dielectroncuts.cfg_max_pt_track < t1.pt()) || (t2.pt() < dielectroncuts.cfg_min_pt_track || dielectroncuts.cfg_max_pt_track < t2.pt())) {
        return false;
      }
      if ((t1.eta() < dielectroncuts.cfg_min_eta_track || dielectroncuts.cfg_max_eta_track < t1.eta()) || (t2.eta() < dielectroncuts.cfg_min_eta_track || dielectroncuts.cfg_max_eta_track < t2.eta())) {
        return false;
      }
      if (v12.Rapidity() < dielectroncuts.cfg_min_pair_y || dielectroncuts.cfg_max_pair_y < v12.Rapidity()) {
        return false;
      }
      return true;
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if ((t1.pt() < dimuoncuts.cfg_min_pt_track || dimuoncuts.cfg_max_pt_track < t1.pt()) || (t2.pt() < dimuoncuts.cfg_min_pt_track || dimuoncuts.cfg_max_pt_track < t2.pt())) {
        return false;
      }
      if ((t1.eta() < dimuoncuts.cfg_min_eta_track || dimuoncuts.cfg_max_eta_track < t1.eta()) || (t2.eta() < dimuoncuts.cfg_min_eta_track || dimuoncuts.cfg_max_eta_track < t2.eta())) {
        return false;
      }
      if (v12.Rapidity() < dimuoncuts.cfg_min_pair_y || dimuoncuts.cfg_max_pair_y < v12.Rapidity()) {
        return false;
      }
      return true;
    } else {
      return false;
    }
    return true;
  }

  template <int sourceId, typename TTrack, typename TMCParticles>
  void fillHistogramsUnfolding(TTrack const& t1, TTrack const& t2, TMCParticles const& mcparticles)
  {
    auto t1mc = mcparticles.iteratorAt(t1.emmcparticleId());
    auto t2mc = mcparticles.iteratorAt(t2.emmcparticleId());

    ROOT::Math::PtEtaPhiMVector v1rec(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2rec(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12rec = v1rec + v2rec;

    ROOT::Math::PtEtaPhiMVector v1mc(t1mc.pt(), t1mc.eta(), t1mc.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2mc(t2mc.pt(), t2mc.eta(), t2mc.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12mc = v1mc + v2mc;

    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[std::make_pair(t1.globalIndex(), t2.globalIndex())];
    }

    if (isPairInAcc(t1, t2) && isPairInAcc(t1mc, t2mc)) { // both rec and mc info are in acceptance.
      if (t1.sign() * t2.sign() < 0) {                    // ULS
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("uls/hsRM"), v12mc.M(), v12mc.Pt(), v12rec.M(), v12rec.Pt(), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("lspp/hsRM"), v12mc.M(), v12mc.Pt(), v12rec.M(), v12rec.Pt(), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("lsmm/hsRM"), v12mc.M(), v12mc.Pt(), v12rec.M(), v12rec.Pt(), weight);
      }
    } else if (!isPairInAcc(t1, t2) && isPairInAcc(t1mc, t2mc)) {
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("uls/hMiss"), v12mc.M(), v12mc.Pt(), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("lspp/hMiss"), v12mc.M(), v12mc.Pt(), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("lsmm/hMiss"), v12mc.M(), v12mc.Pt(), weight);
      }
    } else if (isPairInAcc(t1, t2) && !isPairInAcc(t1mc, t2mc)) {
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("uls/hFake"), v12rec.M(), v12rec.Pt(), weight);
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("lspp/hFake"), v12rec.M(), v12rec.Pt(), weight);
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Unfold/") + HIST(unfolding_dilepton_source_types[sourceId]) + HIST("lsmm/hFake"), v12rec.M(), v12rec.Pt(), weight);
      }
    }
  }

  template <typename TTrack, typename TTracks, typename TCut, typename TMCCollisions, typename TMCParticles>
  bool fillPairUnfolding(TTrack const& t1, TTrack const& t2, TTracks const& tracks, TCut const& cut, TMCCollisions const&, TMCParticles const& mcparticles)
  {
    auto t1mc = mcparticles.iteratorAt(t1.emmcparticleId());
    auto mccollision_from_t1 = t1mc.template emmcevent_as<TMCCollisions>();
    if (cfgEventGeneratorType >= 0 && mccollision_from_t1.getSubGeneratorId() != cfgEventGeneratorType) {
      return false;
    }

    auto t2mc = mcparticles.iteratorAt(t2.emmcparticleId());
    auto mccollision_from_t2 = t2mc.template emmcevent_as<TMCCollisions>();
    if (cfgEventGeneratorType >= 0 && mccollision_from_t2.getSubGeneratorId() != cfgEventGeneratorType) {
      return false;
    }

    if ((std::abs(t1mc.pdgCode()) != pdg_lepton || std::abs(t2mc.pdgCode()) != pdg_lepton) || (t1mc.emmceventId() != t2mc.emmceventId())) {
      return false;
    }
    if (t1mc.pdgCode() * t2mc.pdgCode() > 0) { // ULS
      return false;
    }
    if (!((t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && (t2mc.isPhysicalPrimary() || t2mc.producedByGenerator()))) {
      return false;
    }
    int mother_id = std::max({FindSMULS(t1mc, t2mc, mcparticles), FindSMULS(t2mc, t1mc, mcparticles), FindSMLSPP(t1mc, t2mc, mcparticles), FindSMLSMM(t1mc, t2mc, mcparticles)});
    int hfee_type = IsHF(t1mc, t2mc, mcparticles);
    if (mother_id < 0 && hfee_type < 0) {
      return false;
    }

    if (!isPairOK<true>(t1, t2, cut, tracks)) { // without acceptance
      return false;
    }

    // ROOT::Math::PtEtaPhiMVector v1rec(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    // ROOT::Math::PtEtaPhiMVector v2rec(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    // ROOT::Math::PtEtaPhiMVector v12rec = v1rec + v2rec;

    // ROOT::Math::PtEtaPhiMVector v1mc(t1mc.pt(), t1mc.eta(), t1mc.phi(), leptonM1);
    // ROOT::Math::PtEtaPhiMVector v2mc(t2mc.pt(), t2mc.eta(), t2mc.phi(), leptonM2);
    // ROOT::Math::PtEtaPhiMVector v12mc = v1mc + v2mc;

    if (mother_id > -1) {
      auto mcmother = mcparticles.iteratorAt(mother_id);
      if (!(mcmother.isPhysicalPrimary() || mcmother.producedByGenerator())) {
        return false;
      }
      switch (std::abs(mcmother.pdgCode())) {
        case 111:
        case 221:
        case 331:
        case 113:
        case 223:
        case 333:
        case 443:
        case 100443:
          fillHistogramsUnfolding<0>(t1, t2, mcparticles);
          break;
        default:
          break;
      }
    } else if (hfee_type > -1) {
      switch (hfee_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce):
          fillHistogramsUnfolding<1>(t1, t2, mcparticles);
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be):
          fillHistogramsUnfolding<2>(t1, t2, mcparticles);
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe):
          fillHistogramsUnfolding<2>(t1, t2, mcparticles);
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
          fillHistogramsUnfolding<2>(t1, t2, mcparticles);
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
          fillHistogramsUnfolding<2>(t1, t2, mcparticles);
          break;
        default:
          break;
      }
    }
    return true;
  }

  template <typename TCollisions, typename TTracks1, typename TTracks2, typename TPresilce, typename TCut, typename TAllTracks, typename TMCCollisions, typename TMCParticles>
  void fillUnfolding(TCollisions const& collisions, TTracks1 const& posTracks, TTracks2 const& negTracks, TPresilce const& perCollision, TCut const& cut, TAllTracks const& tracks, TMCCollisions const& mcCollisions, TMCParticles const& mcparticles)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }

      auto posTracks_per_coll = posTracks.sliceByCached(perCollision, collision.globalIndex(), cache); // reconstructed pos tracks
      auto negTracks_per_coll = negTracks.sliceByCached(perCollision, collision.globalIndex(), cache); // reconstructed neg tracks

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillPairUnfolding(pos, neg, tracks, cut, mcCollisions, mcparticles);
      } // end of ULS pairing
      for (const auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        fillPairUnfolding(pos1, pos2, tracks, cut, mcCollisions, mcparticles);
      } // end of LS++ pairing
      for (const auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        fillPairUnfolding(neg1, neg2, tracks, cut, mcCollisions, mcparticles);
      } // end of LS-- pairing
    } // end of collision loop
  }

  SliceCache cache;
  Preslice<MyMCElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Filter trackFilter_electron = nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
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

  Preslice<MyMCMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type;
  Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);

  Filter collisionFilter_centrality = (eventcuts.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventcuts.cfgCentMax) || (eventcuts.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventcuts.cfgCentMax);
  Filter collisionFilter_numContrib = eventcuts.cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < eventcuts.cfgNumContribMax;
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Partition<FilteredMyMCElectrons> positive_electrons = o2::aod::emprimaryelectron::sign > int8_t(0); // reconstructed tracks
  Partition<FilteredMyMCElectrons> negative_electrons = o2::aod::emprimaryelectron::sign < int8_t(0); // reconstructed tracks
  Partition<FilteredMyMCMuons> positive_muons = o2::aod::emprimarymuon::sign > int8_t(0);             // reconstructed tracks
  Partition<FilteredMyMCMuons> negative_muons = o2::aod::emprimarymuon::sign < int8_t(0);             // reconstructed tracks

  Partition<aod::EMMCParticles> positive_electronsMC = o2::aod::mcparticle::pdgCode == -11; // e+
  Partition<aod::EMMCParticles> negative_electronsMC = o2::aod::mcparticle::pdgCode == 11;  // e-
  Partition<aod::EMMCParticles> positive_muonsMC = o2::aod::mcparticle::pdgCode == -13;     // mu+
  Partition<aod::EMMCParticles> negative_muonsMC = o2::aod::mcparticle::pdgCode == 13;      // mu-
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  PresliceUnsorted<aod::EMMCGenVectorMesons> perMcCollision_vm = aod::emmcgenvectormeson::emmceventId;
  // PresliceUnsorted<MyCollisions> recColperMcCollision = aod::emmceventlabel::emmceventId;

  void processAnalysis(FilteredMyCollisions const& collisions, MyMCCollisions const& mccollisions, aod::EMMCParticles const& mcparticles, TLeptons const& leptons)
  {
    // LOGF(info, "collisions.size() = %d, mccollisions.size() = %d, mcparticles.size() = %d", collisions.size(), mccollisions.size(), mcparticles.size());
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, leptons, mccollisions, mcparticles);
      }
      runTruePairing<false>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, leptons, mccollisions, mcparticles);
      runGenInfo<false>(collisions, mccollisions, positive_electronsMC, negative_electronsMC, mcparticles);
      if (cfgFillUnfolding) {
        fillUnfolding(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, leptons, mccollisions, mcparticles);
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, leptons, mccollisions, mcparticles);
      }
      runTruePairing<false>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, leptons, mccollisions, mcparticles);
      runGenInfo<false>(collisions, mccollisions, positive_muonsMC, negative_muonsMC, mcparticles);
      if (cfgFillUnfolding) {
        fillUnfolding(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, leptons, mccollisions, mcparticles);
      }
    }
    map_weight.clear();
  }
  PROCESS_SWITCH(DileptonMC, processAnalysis, "run dilepton mc analysis", true);

  Partition<MySmearedElectrons> positive_electronsMC_smeared = o2::aod::mcparticle::pdgCode == -11; // e+
  Partition<MySmearedElectrons> negative_electronsMC_smeared = o2::aod::mcparticle::pdgCode == 11;  // e-
  Partition<MySmearedMuons> positive_muonsMC_smeared = o2::aod::mcparticle::pdgCode == -13;         // mu+
  Partition<MySmearedMuons> negative_muonsMC_smeared = o2::aod::mcparticle::pdgCode == 13;          // mu-

  void processAnalysis_Smeared(FilteredMyCollisions const& collisions, MyMCCollisions const& mccollisions, TLeptons const& leptons, TSmeardMCParitlces const& mcparticles_smeared)
  {
    // LOGF(info, "collisions.size() = %d, mccollisions.size() = %d, mcparticles.size() = %d", collisions.size(), mccollisions.size(), mcparticles.size());
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, leptons, mccollisions, mcparticles_smeared);
      }
      runTruePairing<true>(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, leptons, mccollisions, mcparticles_smeared);
      runGenInfo<true>(collisions, mccollisions, positive_electronsMC_smeared, negative_electronsMC_smeared, mcparticles_smeared);
      if (cfgFillUnfolding) {
        fillUnfolding(collisions, positive_electrons, negative_electrons, o2::aod::emprimaryelectron::emeventId, fDielectronCut, leptons, mccollisions, mcparticles_smeared);
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (cfgApplyWeightTTCA) {
        fillPairWeightMap(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, leptons, mccollisions, mcparticles_smeared);
      }
      runTruePairing<true>(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, leptons, mccollisions, mcparticles_smeared);
      runGenInfo<true>(collisions, mccollisions, positive_muonsMC_smeared, negative_muonsMC_smeared, mcparticles_smeared);
      if (cfgFillUnfolding) {
        fillUnfolding(collisions, positive_muons, negative_muons, o2::aod::emprimarymuon::emeventId, fDimuonCut, leptons, mccollisions, mcparticles_smeared);
      }
    }
    map_weight.clear();
  }
  PROCESS_SWITCH(DileptonMC, processAnalysis_Smeared, "run dilepton mc analysis with smearing", false);

  void processGen_VM(FilteredMyCollisions const& collisions, MyMCCollisions const&, aod::EMMCGenVectorMesons const& mcparticles)
  {
    // for oemga, phi efficiency
    for (const auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < eventcuts.cfgCentMin || eventcuts.cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      auto mccollision = collision.template emmcevent_as<MyMCCollisions>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }
      auto mctracks_per_coll = mcparticles.sliceBy(perMcCollision_vm, mccollision.globalIndex());

      for (const auto& mctrack : mctracks_per_coll) {
        if (!(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          continue;
        }

        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          if (mctrack.y() < dielectroncuts.cfg_min_pair_y || dielectroncuts.cfg_max_pair_y < mctrack.y()) {
            continue;
          }
        } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
          if (mctrack.y() < dimuoncuts.cfg_min_pair_y || dimuoncuts.cfg_max_pair_y < mctrack.y()) {
            continue;
          }
        }

        switch (std::abs(mctrack.pdgCode())) {
          case 223:
            fRegistry.fill(HIST("Generated/VM/Omega/hPtY"), mctrack.y(), mctrack.pt(), 1.f / mctrack.dsf());
            break;
          case 333:
            fRegistry.fill(HIST("Generated/VM/Phi/hPtY"), mctrack.y(), mctrack.pt(), 1.f / mctrack.dsf());
            break;
          default:
            break;
        }

      } // end of mctracks per mccollision
    } // end of collision loop
  }
  PROCESS_SWITCH(DileptonMC, processGen_VM, "process generated info for vector mesons", false);

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
  PROCESS_SWITCH(DileptonMC, processNorm, "process normalization info", false);

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
  PROCESS_SWITCH(DileptonMC, processBC, "process BC counter", false);

  void processDummy(FilteredMyCollisions const&) {}
  PROCESS_SWITCH(DileptonMC, processDummy, "Dummy function", false);
};

#endif // PWGEM_DILEPTON_CORE_DILEPTONMC_H_
