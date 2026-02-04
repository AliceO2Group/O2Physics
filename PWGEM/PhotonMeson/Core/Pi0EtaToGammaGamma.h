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
/// \file Pi0EtaToGammaGamma.h
/// \brief This code loops over photons and makes pairs for neutral mesons analyses.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMA_H_
#define PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMA_H_

#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/NMHistograms.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
// Dilepton headers
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <EMCALBase/Geometry.h>
#include <EMCALBase/GeometryBase.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/AxisAngle.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

template <o2::aod::pwgem::photonmeson::photonpair::PairType pairtype, o2::soa::is_table... Types>
struct Pi0EtaToGammaGamma {
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  o2::framework::Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  o2::framework::Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};

  o2::framework::Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2"};
  o2::framework::Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  o2::framework::Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  o2::framework::Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  o2::framework::Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  o2::framework::Configurable<float> maxY{"maxY", 0.8, "maximum rapidity for reconstructed particles"};
  o2::framework::Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  o2::framework::Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  o2::framework::ConfigurableAxis ConfVtxBins{"ConfVtxBins", {o2::framework::VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  o2::framework::ConfigurableAxis ConfCentBins{"ConfCentBins", {o2::framework::VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  o2::framework::ConfigurableAxis ConfEPBins{"ConfEPBins", {o2::framework::VARIABLE_WIDTH, -o2::constants::math::PIHalf, -o2::constants::math::PIQuarter, 0.0f, +o2::constants::math::PIQuarter, +o2::constants::math::PIHalf}, "Mixing bins - event plane angle"};
  o2::framework::ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {o2::framework::VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  EMPhotonEventCut fEMEventCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "eventcut_group";
    o2::framework::Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    o2::framework::Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    o2::framework::Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    o2::framework::Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    o2::framework::Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    o2::framework::Configurable<bool> cfgRequireEMCReadoutInMB{"cfgRequireEMCReadoutInMB", false, "require the EMC to be read out in an MB collision (kTVXinEMC)"};
    o2::framework::Configurable<bool> cfgRequireEMCHardwareTriggered{"cfgRequireEMCHardwareTriggered", false, "require the EMC to be hardware triggered (kEMC7 or kDMC7)"};
    o2::framework::Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    o2::framework::Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    o2::framework::Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    o2::framework::Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    o2::framework::Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    o2::framework::Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    o2::framework::Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    o2::framework::Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    o2::framework::Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_max_pt_v0{"cfg_max_pt_v0", 1e+10, "max pT for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.8, "min eta for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.8, "max eta for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    o2::framework::Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    o2::framework::Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    o2::framework::Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    o2::framework::Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min V0 CosPA"};
    o2::framework::Configurable<float> cfg_max_pca{"cfg_max_pca", 1.5, "max distance btween 2 legs"};
    o2::framework::Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    o2::framework::Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    o2::framework::Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to V0"};

    o2::framework::Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    o2::framework::Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    o2::framework::Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    o2::framework::Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};

    o2::framework::Configurable<bool> cfg_apply_ml_cuts{"cfg_apply_ml", false, "flag to apply ML cut"};
    o2::framework::Configurable<bool> cfg_use_2d_binning{"cfg_use_2d_binning", false, "flag to use 2D binning (pT, cent)"};
    o2::framework::Configurable<bool> cfg_load_ml_models_from_ccdb{"cfg_load_ml_models_from_ccdb", true, "flag to load ML models from CCDB"};
    o2::framework::Configurable<int> cfg_timestamp_ccdb{"cfg_timestamp_ccdb", -1, "timestamp for CCDB"};
    o2::framework::Configurable<int> cfg_nclasses_ml{"cfg_nclasses_ml", static_cast<int>(o2::analysis::em_cuts_ml::NCutScores), "number of classes for ML"};
    o2::framework::Configurable<std::string> cfg_cent_type_ml{"cfg_cent_type_ml", "CentFT0C", "centrality type for 2D ML application: CentFT0C, CentFT0M, or CentFT0A"};
    o2::framework::Configurable<std::vector<int>> cfg_cut_dir_ml{"cfg_cut_dir_ml", std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}, "cut direction for ML"};
    o2::framework::Configurable<std::vector<std::string>> cfg_input_feature_names{"cfg_input_feature_names", std::vector<std::string>{"feature1", "feature2"}, "input feature names for ML models"};
    o2::framework::Configurable<std::vector<std::string>> cfg_model_paths_ccdb{"cfg_model_paths_ccdb", std::vector<std::string>{"path_ccdb/BDT_PCM/"}, "CCDB paths for ML models"};
    o2::framework::Configurable<std::vector<std::string>> cfg_onnx_file_names{"cfg_onnx_file_names", std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}, "ONNX file names for ML models"};
    o2::framework::Configurable<std::vector<std::string>> cfg_labels_bins_ml{"cfg_labels_bins_ml", std::vector<std::string>{"bin 0", "bin 1"}, "Labels for bins"};
    o2::framework::Configurable<std::vector<std::string>> cfg_labels_cut_scores_ml{"cfg_labels_cut_scores_ml", std::vector<std::string>{o2::analysis::em_cuts_ml::labelsCutScore}, "Labels for cut scores"};
    o2::framework::Configurable<std::vector<double>> cfg_bins_pt_ml{"cfg_bins_pt_ml", std::vector<double>{0.0, +1e+10}, "pT bin limits for ML application"};
    o2::framework::Configurable<std::vector<double>> cfg_bins_cent_ml{"cfg_bins_cent_ml", std::vector<double>{o2::analysis::em_cuts_ml::vecBinsCent}, "centrality bins for ML"};
    o2::framework::Configurable<std::vector<double>> cfg_cuts_ml_flat{"cfg_cuts_ml_flat", {0.5}, "Flattened ML cuts: [bin0_score0, bin0_score1, ..., binN_scoreM]"};
  } pcmcuts;

  DalitzEECut fDileptonCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    o2::framework::Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    o2::framework::Configurable<float> cfg_max_mass{"cfg_max_mass", 0.1, "max mass"};
    o2::framework::Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    o2::framework::Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    o2::framework::Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    o2::framework::Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    o2::framework::Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    o2::framework::Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    o2::framework::Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.8, "max eta for single track"};
    o2::framework::Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    o2::framework::Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    o2::framework::Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    o2::framework::Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.05, "max dca XY for single track in cm"};
    o2::framework::Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.05, "max dca Z for single track in cm"};
    o2::framework::Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1.5, "max DCA 3D in sigma"};
    o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    o2::framework::Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to electron"};
    o2::framework::Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to enable ITSsa tracks"};
    o2::framework::Configurable<float> cfg_max_pt_track_ITSsa{"cfg_max_pt_track_ITSsa", 0.15, "max pt for ITSsa tracks"};

    o2::framework::Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTOFif), "pid scheme [kTOFif : 0, kTPConly : 1]"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -0.0, "min. TPC n sigma for pion exclusion"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +0.0, "max. TPC n sigma for pion exclusion"};
    o2::framework::Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    o2::framework::Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
  } dileptoncuts;

  EMCPhotonCut fEMCCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "emccut_group";
    o2::framework::Configurable<std::string> clusterDefinition{"clusterDefinition", "kV3Default", "Clusterizer to be selected, e.g. V3Default"};
    o2::framework::Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};
    o2::framework::Configurable<float> EMC_minTime{"EMC_minTime", -20., "Minimum cluster time for EMCal time cut"};
    o2::framework::Configurable<float> EMC_maxTime{"EMC_maxTime", +25., "Maximum cluster time for EMCal time cut"};
    o2::framework::Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    o2::framework::Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    o2::framework::Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    o2::framework::Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    o2::framework::Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    o2::framework::Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    o2::framework::Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    o2::framework::Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    o2::framework::Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
    o2::framework::Configurable<bool> emcUseTM{"emcUseTM", true, "flag to use EMCal track matching cut or not."};
    o2::framework::Configurable<bool> emcUseSecondaryTM{"emcUseSecondaryTM", false, "flag to use EMCal secondary track matching cut or not. Required for PCM-EMC analyses."};
  } emccuts;

  PHOSPhotonCut fPHOSCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "phoscut_group";
    o2::framework::Configurable<float> cfg_min_Ecluster{"cfg_min_Ecluster", 0.3, "Minimum cluster energy for PHOS in GeV"};
  } phoscuts;

  o2::framework::HistogramRegistry fRegistry{"output", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject, false, false};
  // static constexpr std::string_view event_types[2] = {"before/", "after/"};
  // static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  //---------------------------------------------------------------------------
  // Preslices and partitions
  o2::framework::SliceCache cache;
  o2::framework::PresliceOptional<o2::soa::Filtered<o2::soa::Join<o2::aod::V0PhotonsKF, o2::aod::V0KFEMEventIds, o2::aod::V0PhotonsKFPrefilterBitDerived>>> perCollision_pcm = o2::aod::v0photonkf::emeventId;
  o2::framework::PresliceOptional<o2::soa::Join<o2::aod::EmEmcClusters, o2::aod::EMCEMEventIds>> perCollision_emc = o2::aod::emccluster::emeventId;
  o2::framework::PresliceOptional<o2::soa::Join<o2::aod::PHOSClusters, o2::aod::PHOSEMEventIds>> perCollision_phos = o2::aod::phoscluster::emeventId;
  o2::framework::PresliceOptional<o2::soa::Filtered<o2::soa::Join<o2::aod::EMPrimaryElectronsFromDalitz, o2::aod::EMPrimaryElectronEMEventIds, o2::aod::EMPrimaryElectronsPrefilterBitDerived>>> perCollision_electron = o2::aod::emprimaryelectron::emeventId;

  o2::framework::PresliceOptional<o2::aod::EmEmcMTracks> perEMCClusterMT = o2::aod::trackmatching::emEmcClusterId;
  o2::framework::PresliceOptional<o2::aod::EmEmcMSTracks> perEMCClusterMS = o2::aod::trackmatching::emEmcClusterId;

  o2::framework::Partition<o2::soa::Filtered<o2::soa::Join<o2::aod::EMPrimaryElectronsFromDalitz, o2::aod::EMPrimaryElectronEMEventIds, o2::aod::EMPrimaryElectronsPrefilterBitDerived>>> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt&& nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;
  o2::framework::Partition<o2::soa::Filtered<o2::soa::Join<o2::aod::EMPrimaryElectronsFromDalitz, o2::aod::EMPrimaryElectronEMEventIds, o2::aod::EMPrimaryElectronsPrefilterBitDerived>>> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;

  o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, o2::aod::pwgem::dilepton::utils::EMTrack>* emh1 = nullptr;
  o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, o2::aod::pwgem::dilepton::utils::EMTrack>* emh2 = nullptr;
  //---------------------------------------------------------------------------

  std::vector<int> used_photonIds_per_col;                   // <ndf, trackId>
  std::vector<std::pair<int, int>> used_dileptonIds_per_col; // <ndf, trackId>
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  std::vector<float> zvtx_bin_edges;
  std::vector<float> cent_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;

  o2::ccdb::CcdbApi ccdbApi;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  o2::emcal::Geometry* emcalGeom;

  //---------------------------------------------------------------------------
  // In the following are tags defined which help to select the correct preslice and cuts
  struct PCMTag {

    static auto& perCollision()
    {
      static auto slice{o2::aod::v0photonkf::emeventId};
      return slice;
    }

    template <typename Self>
    static auto& cut(Self& s)
    {
      return s.fV0PhotonCut;
    }

    template <typename Self, o2::soa::is_iterator Photon>
    static bool applyCut(Self& s, Photon const& g)
    {
      return s.fV0PhotonCut.template IsSelected<Photon, o2::aod::V0Legs>(g);
    }
  };

  struct EMCTag {

    static auto& perCollision()
    {
      static auto slice{o2::aod::emccluster::emeventId};
      return slice;
    }

    static auto& perClusterMT()
    {
      static auto slice{o2::aod::trackmatching::emEmcClusterId};
      return slice;
    }

    static auto& perClusterMS()
    {
      static auto slice{o2::aod::trackmatching::emEmcClusterId};
      return slice;
    }

    template <typename Self>
    static auto& cut(Self& s)
    {
      return s.fEMCCut;
    }

    // EMCal version has optional tables for matched tracks (global and secondaries)
    template <typename Self, o2::soa::is_iterator Cluster, typename TMatchedTracks = std::nullptr_t, typename TMatchedSecondaries = std::nullptr_t>
    static bool applyCut(Self& s, Cluster const& c, TMatchedTracks const& emcmatchedtracks = nullptr, TMatchedSecondaries const& secondaries = nullptr)
    {
      return s.fEMCCut.IsSelected(c, emcmatchedtracks, secondaries);
    }
  };

  struct PHOSTag {

    static auto& perCollision()
    {
      static auto slice{o2::aod::phoscluster::emeventId};
      return slice;
    }

    template <typename Self>
    static auto& cut(Self& s)
    {
      return s.fPHOSCut;
    }

    template <typename Self, o2::soa::is_iterator Cluster>
    static bool applyCut(Self& s, Cluster const& c)
    {
      return s.fPHOSCut.IsSelected(c);
    }
  };

  struct DalitzEETag {

    static auto& perCollision()
    {
      static auto slice{o2::aod::emprimaryelectron::emeventId};
      return slice;
    }

    template <typename Self>
    static auto& cut(Self& s)
    {
      return s.fDileptonCut;
    }

    // Dalitz version has two tracks as argument + B_z
    template <typename Self, o2::soa::is_iterator Track1, o2::soa::is_iterator Track2>
    static bool applyCut(Self& s, Track1 const& track1, Track2 const& track2, float bz)
    {
      return s.fDileptonCut.IsSelected(track1, track2, bz);
    }
  };

  void init(o2::framework::InitContext&)
  {
    zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
    zvtx_bin_edges.erase(zvtx_bin_edges.begin());

    cent_bin_edges = std::vector<float>(ConfCentBins.value.begin(), ConfCentBins.value.end());
    cent_bin_edges.erase(cent_bin_edges.begin());

    ep_bin_edges = std::vector<float>(ConfEPBins.value.begin(), ConfEPBins.value.end());
    ep_bin_edges.erase(ep_bin_edges.begin());

    LOGF(info, "cfgOccupancyEstimator = %d", cfgOccupancyEstimator.value);
    occ_bin_edges = std::vector<float>(ConfOccupancyBins.value.begin(), ConfOccupancyBins.value.end());
    occ_bin_edges.erase(occ_bin_edges.begin());

    emh1 = new o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, o2::aod::pwgem::dilepton::utils::EMTrack>(ndepth);
    emh2 = new o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, o2::aod::pwgem::dilepton::utils::EMTrack>(ndepth);

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {
      o2::aod::pwgem::photonmeson::utils::nmhistogram::addNMHistograms(&fRegistry, false, "ee#gamma");
    } else {
      o2::aod::pwgem::photonmeson::utils::nmhistogram::addNMHistograms(&fRegistry, false, "#gamma#gamma");
    }
    DefineEMEventCut();
    DefinePCMCut();
    DefineDileptonCut();
    DefineEMCCut();
    DefinePHOSCut();

    if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) {
      fRegistry.addClone("Pair/same/", "Pair/rotation/");
      emcalGeom = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);
    }
    fRegistry.add("Pair/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", o2::framework::kTH1D, {{10001, -0.5, 10000.5}}, true);

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
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
      mRunNumber = collision.runNumber();
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
    fV0PhotonCut.SetD_Bz(d_bz);
    mRunNumber = collision.runNumber();
  }

  ~Pi0EtaToGammaGamma()
  {
    delete emh1;
    emh1 = 0x0;
    delete emh2;
    emh2 = 0x0;

    used_photonIds_per_col.clear();
    used_photonIds_per_col.shrink_to_fit();
    used_dileptonIds_per_col.clear();
    used_dileptonIds_per_col.shrink_to_fit();
    map_mixed_eventId_to_globalBC.clear();
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(eventcuts.cfgZvtxMin, eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireEMCReadoutInMB(eventcuts.cfgRequireEMCReadoutInMB);
    fEMEventCut.SetRequireEMCHardwareTriggered(eventcuts.cfgRequireEMCHardwareTriggered);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, pcmcuts.cfg_max_pt_v0);
    fV0PhotonCut.SetV0EtaRange(pcmcuts.cfg_min_eta_v0, pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfg_max_chi2kf);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfg_max_frac_shared_clusters_tpc);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfg_disable_itsonly_track);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfg_disable_tpconly_track);
    fV0PhotonCut.SetRequireITSTPC(pcmcuts.cfg_require_v0_with_itstpc);
    fV0PhotonCut.SetRequireITSonly(pcmcuts.cfg_require_v0_with_itsonly);
    fV0PhotonCut.SetRequireTPConly(pcmcuts.cfg_require_v0_with_tpconly);

    // for ML
    fV0PhotonCut.SetApplyMlCuts(pcmcuts.cfg_apply_ml_cuts);
    fV0PhotonCut.SetUse2DBinning(pcmcuts.cfg_use_2d_binning);
    fV0PhotonCut.SetLoadMlModelsFromCCDB(pcmcuts.cfg_load_ml_models_from_ccdb);
    fV0PhotonCut.SetNClassesMl(pcmcuts.cfg_nclasses_ml);
    fV0PhotonCut.SetMlTimestampCCDB(pcmcuts.cfg_timestamp_ccdb);
    fV0PhotonCut.SetCcdbUrl(ccdburl);
    fV0PhotonCut.SetCentralityTypeMl(pcmcuts.cfg_cent_type_ml);
    fV0PhotonCut.SetCutDirMl(pcmcuts.cfg_cut_dir_ml);
    fV0PhotonCut.SetMlModelPathsCCDB(pcmcuts.cfg_model_paths_ccdb);
    fV0PhotonCut.SetMlOnnxFileNames(pcmcuts.cfg_onnx_file_names);
    fV0PhotonCut.SetBinsPtMl(pcmcuts.cfg_bins_pt_ml);
    fV0PhotonCut.SetBinsCentMl(pcmcuts.cfg_bins_cent_ml);
    fV0PhotonCut.SetCutsMl(pcmcuts.cfg_cuts_ml_flat);
    fV0PhotonCut.SetNamesInputFeatures(pcmcuts.cfg_input_feature_names);
    fV0PhotonCut.SetLabelsBinsMl(pcmcuts.cfg_labels_bins_ml);
    fV0PhotonCut.SetLabelsCutScoresMl(pcmcuts.cfg_labels_cut_scores_ml);

    if (pcmcuts.cfg_apply_ml_cuts) {
      fV0PhotonCut.initV0MlModels(ccdbApi);
    }
  }

  void DefineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mass, dileptoncuts.cfg_max_mass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.ApplyPhiV(dileptoncuts.cfg_apply_phiv);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    // for track
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_max_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetMaxFracSharedClustersTPC(dileptoncuts.cfg_max_frac_shared_clusters_tpc);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);
    fDileptonCut.SetTrackDca3DRange(0.f, dileptoncuts.cfg_max_dca3dsigma_track); // in sigma
    fDileptonCut.IncludeITSsa(dileptoncuts.includeITSsa, dileptoncuts.cfg_max_pt_track_ITSsa);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);
  }

  void DefineEMCCut()
  {
    fEMCCut = EMCPhotonCut("fEMCCut", "fEMCCut");

    fEMCCut.SetClusterizer(emccuts.clusterDefinition);
    fEMCCut.SetMinE(emccuts.EMC_minE);
    fEMCCut.SetMinNCell(emccuts.EMC_minNCell);
    fEMCCut.SetM02Range(emccuts.EMC_minM02, emccuts.EMC_maxM02);
    fEMCCut.SetTimeRange(emccuts.EMC_minTime, emccuts.EMC_maxTime);

    fEMCCut.SetTrackMatchingEtaParams(emccuts.EMC_TM_Eta->at(0), emccuts.EMC_TM_Eta->at(1), emccuts.EMC_TM_Eta->at(2));
    fEMCCut.SetTrackMatchingPhiParams(emccuts.EMC_TM_Phi->at(0), emccuts.EMC_TM_Phi->at(1), emccuts.EMC_TM_Phi->at(2));

    fEMCCut.SetMinEoverP(emccuts.EMC_Eoverp);
    fEMCCut.SetUseExoticCut(emccuts.EMC_UseExoticCut);
    fEMCCut.SetUseTM(emccuts.emcUseTM.value);                   // disables or enables TM
    fEMCCut.SetUseSecondaryTM(emccuts.emcUseSecondaryTM.value); // disables or enables secondary TM
  }

  void DefinePHOSCut()
  {
    fPHOSCut.SetEnergyRange(phoscuts.cfg_min_Ecluster, 1e+10);
  }

  /// \brief returns if cluster is too close to edge of EMCal (using rotation background method only for EMCal!)
  bool IsTooCloseToEdge(const int cellID, const int DistanceToBorder = 1)
  {
    if (DistanceToBorder <= 0) {
      return false;
    }
    if (cellID < 0) {
      return true;
    }

    int iBadCell = -1;

    // check distance to border in case the cell is okay
    auto [iSupMod, iMod, iPhi, iEta] = emcalGeom->GetCellIndex(cellID);
    auto [irow, icol] = emcalGeom->GetCellPhiEtaIndexInSModule(iSupMod, iMod, iPhi, iEta);

    // Check rows/phi
    int iRowLast = 24;
    if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_HALF) {
      iRowLast /= 2; // 2/3 sm case
    } else if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::EMCAL_THIRD) {
      iRowLast /= 3; // 1/3 sm case
    } else if (emcalGeom->GetSMType(iSupMod) == o2::emcal::EMCALSMType::DCAL_EXT) {
      iRowLast /= 3; // 1/3 sm case
    }

    if (irow < DistanceToBorder || (iRowLast - irow) <= DistanceToBorder) {
      iBadCell = 1;
    }

    if (iBadCell > 0) {
      return true;
    }
    return false;
  }

  /// \brief Calculate background (using rotation background method only for EMCal!)
  template <typename TPhotons>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, float eventWeight)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < 3) {
      return;
    }
    const float rotationAngle = o2::constants::math::PIHalf; // rotaion angle 90 degree
    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationAngle);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    int iCellID_photon1 = 0;
    int iCellID_photon2 = 0;

    try {
      iCellID_photon1 = emcalGeom->GetAbsCellIdFromEtaPhi(photon1.Eta(), photon1.Phi());
      if (IsTooCloseToEdge(iCellID_photon1, emccuts.cfgDistanceToEdge.value)) {
        iCellID_photon1 = -1;
      }
    } catch (o2::emcal::InvalidPositionException& e) {
      iCellID_photon1 = -1;
    }
    try {
      iCellID_photon2 = emcalGeom->GetAbsCellIdFromEtaPhi(photon2.Eta(), photon2.Phi());
      if (IsTooCloseToEdge(iCellID_photon2, emccuts.cfgDistanceToEdge.value)) {
        iCellID_photon2 = -1;
      }
    } catch (o2::emcal::InvalidPositionException& e) {
      iCellID_photon2 = -1;
    }
    if (iCellID_photon1 == -1 && iCellID_photon2 == -1) {
      return;
    }

    for (const auto& photon : photons_coll) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      if (!(fEMCCut.IsSelected(photon))) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector photon3(photon.pt(), photon.eta(), photon.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
      ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;

      float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));
      float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));

      if (openingAngle1 > emccuts.minOpenAngle && std::fabs(mother1.Rapidity()) < maxY && iCellID_photon1 > 0) {
        fRegistry.fill(HIST("Pair/rotation/hs"), mother1.M(), mother1.Pt(), eventWeight);
      }
      if (openingAngle2 > emccuts.minOpenAngle && std::fabs(mother2.Rapidity()) < maxY && iCellID_photon2 > 0) {
        fRegistry.fill(HIST("Pair/rotation/hs"), mother2.M(), mother2.Pt(), eventWeight);
      }
    }
    return;
  }

  /// \brief function to run the photon pairing
  /// \tparam TDetectorTag1 tag for TPhotons1 type to select the proper cut function and arguments
  /// \tparam TDetectorTag2 tag for TPhotons2 type to select the proper cut function and arguments
  /// \tparam TCombinationPolicy StrictlyUpperIndex or FullIndex depending if we combine from the same detector or different detectors
  /// \tparam TCollisions collision table type
  /// \tparam TPhotons1 first photon table type
  /// \tparam TPhotons2 second photon table type
  /// \tparam TLegs V0 leg table type used in TCut1 and/or TCut2 (only for PCM/DalitzEE)
  /// \tparam TMatchedTracks matched global track table type for EMCal (only for EMC)
  /// \tparam TMatchedSecondaries matched secondary track table type for EMCal (only for EMC)
  /// \param collisions table of collisions (TCollisions)
  /// \param photons1 table of photons (TPhotons1)
  /// \param photons2 table of photons (TPhotons2)
  /// \param legs placeholder argument used only for template deduction (PCM/DalitzEE)
  /// \param matchedtracks table of matched global tracks to EMCal clusters (optional)
  /// \param matchedsecondaries table of matched secondary tracks to EMCal clusters (optional)
  template <typename TDetectorTag1, typename TDetectorTag2, template <typename...> class TCombinationPolicy = o2::soa::CombinationsStrictlyUpperIndexPolicy, o2::soa::is_table TCollisions, o2::soa::is_table TPhotons1, o2::soa::is_table TPhotons2, typename TLegs = std::nullptr_t, typename TMatchedTracks = std::nullptr_t, typename TMatchedSecondaries = std::nullptr_t>
  void runPairing(TCollisions const& collisions,
                  TPhotons1 const& photons1, TPhotons2 const& photons2, TLegs const& /*legs*/ = nullptr, TMatchedTracks const& matchedTracks = nullptr, TMatchedSecondaries const& matchedSecondaries = nullptr)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      fV0PhotonCut.SetCentrality(collision.centFT0A(), collision.centFT0C(), collision.centFT0M());
      int ndiphoton = 0;
      if ((pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS || pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }

      float weight = 1.f;
      if constexpr (std::is_same_v<std::decay_t<TCollisions>, o2::soa::Filtered<o2::soa::Join<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMEventsQvec>, o2::aod::EMEventsWeight>>>) {
        weight = collision.weight();
      }

      if (eventcuts.onlyKeepWeightedEvents && std::fabs(weight - 1.f) < 1E-10) {
        continue;
      }

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, weight);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, weight);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0, weight); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0, weight);  // accepted

      int zbin = lower_bound(zvtx_bin_edges.begin(), zvtx_bin_edges.end(), collision.posZ()) - zvtx_bin_edges.begin() - 1;
      if (zbin < 0) {
        zbin = 0;
      } else if (static_cast<int>(zvtx_bin_edges.size()) - 2 < zbin) {
        zbin = static_cast<int>(zvtx_bin_edges.size()) - 2;
      }

      float centrality = centralities[cfgCentEstimator];
      int centbin = lower_bound(cent_bin_edges.begin(), cent_bin_edges.end(), centrality) - cent_bin_edges.begin() - 1;
      if (centbin < 0) {
        centbin = 0;
      } else if (static_cast<int>(cent_bin_edges.size()) - 2 < centbin) {
        centbin = static_cast<int>(cent_bin_edges.size()) - 2;
      }

      float ep2 = collision.ep2btot();
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
      std::pair<int, int64_t> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {
        auto photons1_per_collision = photons1.sliceByCached(TDetectorTag1::perCollision(), collision.globalIndex(), cache);
        auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

        for (const auto& g1 : photons1_per_collision) {
          if constexpr (std::is_same_v<TDetectorTag1, PCMTag>) {
            if (!TDetectorTag1::applyCut(*this, g1)) {
              continue;
            }
          } else {
            if (!TDetectorTag2::applyCut(*this, g1)) {
              continue;
            }
          }
          auto pos1 = g1.template posTrack_as<TLegs>();
          auto ele1 = g1.template negTrack_as<TLegs>();
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

          for (const auto& [pos2, ele2] : o2::soa::combinations(TCombinationPolicy(positrons_per_collision, electrons_per_collision))) {

            if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks.
              continue;
            }
            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
              continue;
            }

            if constexpr (std::is_same_v<TDetectorTag2, DalitzEETag>) {
              if (!TDetectorTag2::applyCut(*this, pos2, ele2, d_bz)) {
                continue;
              }
            } else {
              if (!TDetectorTag1::applyCut(*this, pos2, ele2, d_bz)) {
                continue;
              }
            }
            // if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
            //   continue;
            // }

            // if (!cut2.IsSelectedPair(pos2, ele2, d_bz)) {
            //   continue;
            // }

            ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
            ROOT::Math::PtEtaPhiMVector veeg = v_gamma + v_pos + v_ele;
            if (std::fabs(veeg.Rapidity()) > maxY) {
              continue;
            }

            fRegistry.fill(HIST("Pair/same/hs"), veeg.M(), veeg.Pt(), weight);

            std::pair<int, int> tuple_tmp_id2 = std::make_pair(pos2.trackId(), ele2.trackId());
            if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g1.globalIndex()) == used_photonIds_per_col.end()) {
              emh1->AddTrackToEventPool(key_df_collision, o2::aod::pwgem::dilepton::utils::EMTrack(g1.pt(), g1.eta(), g1.phi(), 0));
              used_photonIds_per_col.emplace_back(g1.globalIndex());
            }
            if (std::find(used_dileptonIds_per_col.begin(), used_dileptonIds_per_col.end(), tuple_tmp_id2) == used_dileptonIds_per_col.end()) {
              emh2->AddTrackToEventPool(key_df_collision, o2::aod::pwgem::dilepton::utils::EMTrack(v_ee.Pt(), v_ee.Eta(), v_ee.Phi(), v_ee.M()));
              used_dileptonIds_per_col.emplace_back(tuple_tmp_id2);
            }
            ndiphoton++;
          } // end of dielectron loop
        } // end of g1 loop
      } else { // PCM-PCM, EMC-EMC, PHOS-PHOS, PCM-EMC and PCM-PHOS.
        auto photons1_per_collision = photons1.sliceByCached(TDetectorTag1::perCollision(), collision.globalIndex(), cache);
        auto photons2_per_collision = photons2.sliceByCached(TDetectorTag2::perCollision(), collision.globalIndex(), cache);

        for (const auto& [g1, g2] : o2::soa::combinations(TCombinationPolicy(photons1_per_collision, photons2_per_collision))) {
          if constexpr (std::is_same_v<TDetectorTag1, EMCTag>) {
            // For the EMCal case we need to get the primary and secondary matched tracks
            auto matchedTracks1 = matchedTracks.sliceByCached(TDetectorTag1::perClusterMT(), g1.globalIndex(), cache);
            auto matchedSecondaries1 = matchedSecondaries.sliceByCached(TDetectorTag1::perClusterMS(), g1.globalIndex(), cache);

            if (!TDetectorTag1::applyCut(*this, g1, matchedTracks1, matchedSecondaries1)) {
              continue;
            }
          } else {
            if (!TDetectorTag1::applyCut(*this, g1) || !TDetectorTag2::applyCut(*this, g2)) {
              continue;
            }
          }

          if constexpr (std::is_same_v<TDetectorTag2, EMCTag>) {
            auto matchedTracks2 = matchedTracks.sliceByCached(TDetectorTag2::perClusterMT(), g2.globalIndex(), cache);
            auto matchedSecondaries2 = matchedSecondaries.sliceByCached(TDetectorTag2::perClusterMS(), g2.globalIndex(), cache);
            if (!TDetectorTag2::applyCut(*this, g2, matchedTracks2, matchedSecondaries2)) {
              continue;
            }
          } else {
            if (!TDetectorTag2::applyCut(*this, g2)) {
              continue;
            }
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (std::fabs(v12.Rapidity()) > maxY) {
            continue;
          }

          fRegistry.fill(HIST("Pair/same/hs"), v12.M(), v12.Pt(), weight);

          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g1.globalIndex()) == used_photonIds_per_col.end()) {
            emh1->AddTrackToEventPool(key_df_collision, o2::aod::pwgem::dilepton::utils::EMTrack(g1.pt(), g1.eta(), g1.phi(), 0));
            used_photonIds_per_col.emplace_back(g1.globalIndex());
          }
          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g2.globalIndex()) == used_photonIds_per_col.end()) {
            emh2->AddTrackToEventPool(key_df_collision, o2::aod::pwgem::dilepton::utils::EMTrack(g2.pt(), g2.eta(), g2.phi(), 0));
            used_photonIds_per_col.emplace_back(g2.globalIndex());
          }
          ndiphoton++;
        } // end of pairing loop
      } // end of pairing in same event

      used_photonIds_per_col.clear();
      used_photonIds_per_col.shrink_to_fit();
      used_dileptonIds_per_col.clear();
      used_dileptonIds_per_col.shrink_to_fit();

      // event mixing
      if (!cfgDoMix || !(ndiphoton > 0)) {
        continue;
      }

      // make a vector of selected photons in this collision.
      auto selected_photons1_in_this_event = emh1->GetTracksPerCollision(key_df_collision);
      auto selected_photons2_in_this_event = emh2->GetTracksPerCollision(key_df_collision);

      auto collisionIds1_in_mixing_pool = emh1->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIds2_in_mixing_pool = emh2->GetCollisionIdsFromEventPool(key_bin);

      if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM || pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS || pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) { // same kinds pairing
        for (const auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), ngamma = %d | event pool (%d, %d), ngamma = %d", ndf, collision.globalIndex(), selected_photons1_in_this_event.size(), mix_dfId, mix_collisionId, photons1_from_event_pool.size());

          for (const auto& g1 : selected_photons1_in_this_event) {
            for (const auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (std::fabs(v12.Rapidity()) > maxY) {
                continue;
              }

              fRegistry.fill(HIST("Pair/mix/hs"), v12.M(), v12.Pt(), weight);
            }
          }
        } // end of loop over mixed event pool

      } else { // [photon1 from event1, photon2 from event2] and [photon1 from event2, photon2 from event1]
        for (const auto& mix_dfId_collisionId : collisionIds2_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto photons2_from_event_pool = emh2->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), ngamma = %d | event pool (%d, %d), nll = %d", ndf, collision.globalIndex(), selected_photons1_in_this_event.size(), mix_dfId, mix_collisionId, photons2_from_event_pool.size());

          for (const auto& g1 : selected_photons1_in_this_event) {
            for (const auto& g2 : photons2_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v2.SetM(g2.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (std::fabs(v12.Rapidity()) > maxY) {
                continue;
              }
              fRegistry.fill(HIST("Pair/mix/hs"), v12.M(), v12.Pt(), weight);
            }
          }
        } // end of loop over mixed event pool
        for (const auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), nll = %d | event pool (%d, %d), ngamma = %d", ndf, collision.globalIndex(), selected_photons2_in_this_event.size(), mix_dfId, mix_collisionId, photons1_from_event_pool.size());

          for (const auto& g1 : selected_photons2_in_this_event) {
            for (const auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v1.SetM(g1.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (std::fabs(v12.Rapidity()) > maxY) {
                continue;
              }
              fRegistry.fill(HIST("Pair/mix/hs"), v12.M(), v12.Pt(), weight);
            }
          }
        } // end of loop over mixed event pool
      }

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh2->AddCollisionIdAtLast(key_bin, key_df_collision);
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
      }

    } // end of collision loop
  }

  o2::framework::expressions::Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  o2::framework::expressions::Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  o2::framework::expressions::Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  // using FilteredMyCollisions = o2::soa::Filtered<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMEventsQvec>>;

  o2::framework::expressions::Filter prefilter_pcm = ifnode(pcmcuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::v0photonkf::pfbderived == static_cast<uint16_t>(0), true);
  o2::framework::expressions::Filter prefilter_primaryelectron = ifnode(dileptoncuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::emprimaryelectron::pfbderived == static_cast<uint16_t>(0), true);

  int ndf = 0;
  void processAnalysis(o2::soa::Filtered<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMEventsQvec>> const& collisions, Types const&... args)
  {
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM) {
      auto&& [v0photons, v0legs] = std::forward_as_tuple(args...);
      // auto v0photons = std::get<0>(std::tie(args...));
      // auto v0legs = std::get<1>(std::tie(args...));
      runPairing<PCMTag, PCMTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, v0photons, v0photons, v0legs);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {
      auto&& [v0photons, v0legs, emprimaryelectrons] = std::forward_as_tuple(args...);
      // auto v0photons = std::get<0>(std::tie(args...));
      // auto v0legs = std::get<1>(std::tie(args...));
      // auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runPairing<PCMTag, DalitzEETag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, emprimaryelectrons, v0legs);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) {
      auto&& [emcClusters, emcMatchedTracks, emcMatchedSecondaries] = std::forward_as_tuple(args...);
      // auto emcClusters = std::get<0>(std::tie(args...));
      // auto emcMatchedTracks = std::get<1>(std::tie(args...));
      // auto emcMatchedSecondaries = std::get<2>(std::tie(args...));
      runPairing<EMCTag, EMCTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, emcClusters, emcClusters, nullptr, emcMatchedTracks, emcMatchedSecondaries);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS) {
      auto&& [phosclusters] = std::forward_as_tuple(args...);
      // auto phosclusters = std::get<0>(std::tie(args...));
      runPairing<PHOSTag, PHOSTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, phosclusters, phosclusters);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMEMC) {
      auto&& [v0photons, emcclusters, v0legs, emcMatchedTracks, emcMatchedSecondaries] = std::forward_as_tuple(args...);
      runPairing<PCMTag, EMCTag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, emcclusters, v0legs, emcMatchedTracks, emcMatchedSecondaries);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPHOS) {
      auto&& [v0photons, v0legs, phosclusters] = std::forward_as_tuple(args...);
      runPairing<PCMTag, PHOSTag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, phosclusters, v0legs);
    }
    ndf++;
  }
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processAnalysis, "process pair analysis", true);

  // using FilteredMyCollisionsWithJJMC = o2::soa::Filtered<o2::soa::Join<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMEventsQvec>, o2::aod::EMEventsWeight>>;
  void processAnalysisJJMC(o2::soa::Filtered<o2::soa::Join<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMEventsQvec>, o2::aod::EMEventsWeight>> const& collisions, Types const&... args)
  {
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM) {
      auto&& [v0photons, v0legs] = std::forward_as_tuple(args...);
      runPairing<PCMTag, PCMTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, v0photons, v0photons, v0legs);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {
      auto&& [v0photons, v0legs, emprimaryelectrons] = std::forward_as_tuple(args...);
      runPairing<PCMTag, DalitzEETag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, emprimaryelectrons, v0legs);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) {
      auto&& [emcClusters, emcMatchedTracks, emcMatchedSecondaries] = std::forward_as_tuple(args...);
      runPairing<EMCTag, EMCTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, emcClusters, emcClusters, nullptr, emcMatchedTracks, emcMatchedSecondaries);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS) {
      auto&& [phosclusters] = std::forward_as_tuple(args...);
      // auto phosclusters = std::get<0>(std::tie(args...));
      runPairing<PHOSTag, PHOSTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, phosclusters, phosclusters);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMEMC) {
      auto&& [v0photons, emcclusters, v0legs, emcMatchedTracks, emcMatchedSecondaries] = std::forward_as_tuple(args...);
      runPairing<PCMTag, EMCTag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, emcclusters, v0legs, emcMatchedTracks, emcMatchedSecondaries);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPHOS) {
      auto&& [v0photons, phosclusters, v0legs] = std::forward_as_tuple(args...);
      runPairing<PCMTag, PHOSTag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, phosclusters, v0legs);
    }
    ndf++;
  }
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processAnalysisJJMC, "process pair analysis", false);

  void processDummy(o2::aod::EMEvents const&) {}
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processDummy, "Dummy function", false);
};
#endif // PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMA_H_
