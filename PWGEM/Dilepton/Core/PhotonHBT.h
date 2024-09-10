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
// This code loops over v0 photons and makes pairs for photon HBT analysis.
//    Please write to: daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_CORE_PHOTONHBT_H_
#define PWGEM_DILEPTON_CORE_PHOTONHBT_H_

#include <algorithm>
#include <cstring>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <tuple>
#include <vector>

#include "TString.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"

namespace o2::aod::pwgem::dilepton::core::photonhbt
{
enum class ggHBTPairType : int {
  kPCMPCM = 0,
  kPCMEE = 1,
  kEEEE = 2,
};
} // namespace o2::aod::pwgem::dilepton::core::photonhbt

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils;
using namespace o2::aod::pwgem::dilepton::core::photonhbt;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithSWT = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMSWTriggerInfos>;
using MyCollisionWithSWT = MyCollisionsWithSWT::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;
using FilteredMyTracks = soa::Filtered<MyTracks>;
using FilteredMyTrack = FilteredMyTracks::iterator;

template <ggHBTPairType pairtype, typename... Types>
struct PhotonHBT {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<bool> cfgDo3D{"cfgDo3D", false, "enable 3D analysis"};
  Configurable<int> cfgEP2Estimator_for_Mix{"cfgEP2Estimator_for_Mix", 3, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2, NTPV:3"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  // Configurable<float> cfgSpherocityMin{"cfgSpherocityMin", -999.f, "min. spherocity"};
  // Configurable<float> cfgSpherocityMax{"cfgSpherocityMax", +999.f, "max. spherocity"};
  Configurable<float> maxY{"maxY", 0.8, "maximum rapidity for reconstructed particles"};
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
  Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", true, "flag to apply weighting by 1/N"};

  ConfigurableAxis ConfKtBins{"ConfKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}, "kT bins for output histograms"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -1, "min. occupancy"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<bool> cfg_require_v0_on_wwire_ib{"cfg_require_v0_on_wwire_ib", false, "flag to select V0s on W wires ITSib"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.997, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
    Configurable<bool> cfg_require_v0_with_correct_xz{"cfg_require_v0_with_correct_xz", true, "flag to select V0s with correct xz"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};

    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 10, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
  } pcmcuts;

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "dielectroncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.015, "max mass"}; // this is valid, because only ULS is used.
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pT"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pT"};
    Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -0.8, "min pair rapidity"};
    Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", +0.8, "max pair rapidity"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 2.0, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "max eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<float> cfg_min_its_cluster_size{"cfg_min_its_cluster_size", 0.f, "min ITS cluster size"};
    Configurable<float> cfg_max_its_cluster_size{"cfg_max_its_cluster_size", 16.f, "max ITS cluster size"};
    Configurable<float> cfg_max_p_its_cluster_size{"cfg_max_p_its_cluster_size", 0.2, "max p to apply ITS cluster size cut"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -3.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_pin_pirejTPC{"cfg_max_pin_pirejTPC", 0.5, "max. pin for pion rejection in TPC"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};

    // CCDB configuration for PID ML
    Configurable<std::string> BDTLocalPathGamma{"BDTLocalPathGamma", "pid_ml_xgboost.onnx", "Path to the local .onnx file"};

    Configurable<std::string> BDTPathCCDB{"BDTPathCCDB", "Users/d/dsekihat/pwgem/pidml/", "Path on CCDB"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dielectroncuts;

  struct : ConfigurableGroup {
    std::string prefix = "ggpaircut_group";
    Configurable<bool> applydR{"applydR", false, "apply deta-dphi cut to avoid track splitting/merging"};
    Configurable<float> cfgMinDeltaEta{"cfgMinDeltaEta", 0.f, "min. delta-eta between 2 photons"};
    Configurable<float> cfgMinDeltaPhi{"cfgMinDeltaPhi", 0.f, "min. delta-phi between 2 photons"};
    Configurable<float> cfgMinDeltaR{"cfgMinDeltaR", 0.f, "min. delta-r between 2 photons"};
    Configurable<float> cfgMinDeltaZ{"cfgMinDeltaZ", 0.f, "min. delta-z between 2 photons"};
  } ggpaircuts;

  ~PhotonHBT()
  {
    delete emh1;
    emh1 = 0x0;
    delete emh2;
    emh2 = 0x0;

    map_mixed_eventId_to_globalBC.clear();

    used_photonIds.clear();
    used_photonIds.shrink_to_fit();
    used_dileptonIds.clear();
    used_dileptonIds.shrink_to_fit();

    if (eid_bdt) {
      delete eid_bdt;
    }
  }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before", "after"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  std::vector<float> zvtx_bin_edges;
  std::vector<float> cent_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;

  void init(InitContext& /*context*/)
  {
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

    emh1 = new MyEMH(ndepth);
    emh2 = new MyEMH(ndepth);

    DefineEMEventCut();
    DefinePCMCut();
    DefineDileptonCut();

    addhistograms();
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

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

    if constexpr (isTriggerAnalysis) {
      LOGF(info, "Trigger analysis is enabled. Desired trigger name = %s", cfg_swt_name.value);
      LOGF(info, "total inspected TVX events = %d in run number %d", collision.nInspectedTVX(), collision.runNumber());
      fRegistry.fill(HIST("Event/hNInspectedTVX"), collision.runNumber(), collision.nInspectedTVX());
    }
  }

  void addhistograms()
  {
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);
    std::string_view qvec_det_names[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix", Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEP2Estimator_for_Mix].data()), kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry.add("Event/after/hEP2_CentFT0C_forMix", Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEP2Estimator_for_Mix].data()), kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);

    // pair info
    const AxisSpec axis_kt{ConfKtBins, "k_{T} (GeV/c)"};
    const AxisSpec axis_qinv{60, 0.0, +0.3, "q_{inv} (GeV/c)"};
    const AxisSpec axis_kstar{60, 0.0, +0.3, "k* (GeV/c)"};
    const AxisSpec axis_qabs_lcms{60, 0.0, +0.3, "|#bf{q}|^{LCMS} (GeV/c)"};
    const AxisSpec axis_qout{60, 0.0, +0.3, "q_{out} (GeV/c)"};   // qout does not change between LAB and LCMS frame
    const AxisSpec axis_qside{60, 0.0, +0.3, "q_{side} (GeV/c)"}; // qside does not change between LAB and LCMS frame
    const AxisSpec axis_qlong{60, 0.0, +0.3, "q_{long} (GeV/c)"};

    if (cfgDo3D) {
      fRegistry.add("Pair/same/hs_3d", "diphoton correlation 3D LCMS", kTHnSparseD, {axis_qout, axis_qside, axis_qlong, axis_kt}, true);
    } else {
      if constexpr (pairtype == ggHBTPairType::kPCMPCM) { // identical particle femtoscopy
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D", kTHnSparseD, {axis_qinv, axis_qabs_lcms, axis_kt}, true);
      } else { // non-identical particle femtoscopy
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D", kTHnSparseD, {axis_kstar, axis_qabs_lcms, axis_kt}, true);
      }
    }

    if constexpr (pairtype == ggHBTPairType::kPCMPCM) { // dr, dz of conversion points
      fRegistry.add("Pair/same/hDeltaRDeltaZ", "diphoton distance in RZ;#Deltar = #sqrt{(#Deltax)^{2} + (#Deltay)^{2}} (cm);#Deltaz (cm)", kTH2D, {{100, 0, +10}, {200, -10, 10}}, true);
    } else { // deta, dphi of track momentum
      fRegistry.add("Pair/same/hDeltaEtaDeltaPhi", "diphoton distance in #eta-#varphi plane;#Delta#varphi (rad.);#Delta#eta", kTH2D, {{200, -0.1, +0.1}, {200, -0.1, 0.1}}, true);
    }

    fRegistry.addClone("Pair/same/", "Pair/mix/");
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

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, 1e10f);
    fV0PhotonCut.SetV0EtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetTrackPtRange(pcmcuts.cfg_min_pt_v0 * 0.5, 1e+10f);
    fV0PhotonCut.SetTrackEtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfg_disable_itsonly_track);

    if (pcmcuts.cfg_reject_v0_on_itsib) {
      fV0PhotonCut.SetNClustersITS(2, 4);
    } else {
      fV0PhotonCut.SetNClustersITS(0, 7);
    }
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetIsWithinBeamPipe(pcmcuts.cfg_require_v0_with_correct_xz);

    if (pcmcuts.cfg_require_v0_with_itstpc) {
      fV0PhotonCut.SetRequireITSTPC(true);
      fV0PhotonCut.SetMaxPCA(1.0);
      fV0PhotonCut.SetRxyRange(4, 40);
    }
    if (pcmcuts.cfg_require_v0_with_itsonly) {
      fV0PhotonCut.SetRequireITSonly(true);
      fV0PhotonCut.SetMaxPCA(1.0);
      fV0PhotonCut.SetRxyRange(4, 24);
    }
    if (pcmcuts.cfg_require_v0_with_tpconly) {
      fV0PhotonCut.SetRequireTPConly(true);
      fV0PhotonCut.SetMaxPCA(3.0);
      fV0PhotonCut.SetRxyRange(36, 90);
    }
    if (pcmcuts.cfg_require_v0_on_wwire_ib) {
      fV0PhotonCut.SetMaxPCA(0.3);
      fV0PhotonCut.SetOnWwireIB(true);
      fV0PhotonCut.SetOnWwireOB(false);
      fV0PhotonCut.SetRxyRange(7, 14);
    }
  }

  o2::ml::OnnxModel* eid_bdt = nullptr;
  void DefineDileptonCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for pair
    fDielectronCut.SetMeeRange(dielectroncuts.cfg_min_mass, dielectroncuts.cfg_max_mass);
    fDielectronCut.SetPairPtRange(dielectroncuts.cfg_min_pair_pt, dielectroncuts.cfg_max_pair_pt);
    fDielectronCut.SetPairYRange(dielectroncuts.cfg_min_pair_y, dielectroncuts.cfg_max_pair_y);
    fDielectronCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dielectroncuts.cfg_phiv_intercept) / dielectroncuts.cfg_phiv_slope; });
    fDielectronCut.SetPairDCARange(dielectroncuts.cfg_min_pair_dca3d, dielectroncuts.cfg_max_pair_dca3d); // in sigma
    fDielectronCut.ApplyPhiV(dielectroncuts.cfg_apply_phiv);
    fDielectronCut.ApplyPrefilter(dielectroncuts.cfg_apply_pf);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, 1e+10f);
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITS(dielectroncuts.cfg_min_its_cluster_size, dielectroncuts.cfg_max_its_cluster_size, dielectroncuts.cfg_max_p_its_cluster_size);
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
    fDielectronCut.SetMaxPinForPionRejectionTPC(dielectroncuts.cfg_max_pin_pirejTPC);

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

  template <int ev_id, typename TCollision>
  void fillPairHistogram(TCollision const&, const ROOT::Math::PtEtaPhiMVector v1, const ROOT::Math::PtEtaPhiMVector v2, const float weight = 1.f)
  {
    // Lab. frame
    ROOT::Math::PtEtaPhiMVector q12 = v1 - v2;
    ROOT::Math::PtEtaPhiMVector k12 = 0.5 * (v1 + v2);
    float qinv = -q12.M(); // for identical particles -> qinv = 2 x kstar
    float kt = k12.Pt();

    // ROOT::Math::XYZVector q_3d = q12.Vect();                                   // 3D q vector
    ROOT::Math::XYZVector uv_out(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0); // unit vector for out. i.e. parallel to kt
    ROOT::Math::XYZVector uv_long(0, 0, 1);                                    // unit vector for long, beam axis
    ROOT::Math::XYZVector uv_side = uv_out.Cross(uv_long);                     // unit vector for side

    ROOT::Math::PxPyPzEVector v1_cartesian(v1);
    ROOT::Math::PxPyPzEVector v2_cartesian(v2);
    ROOT::Math::PxPyPzEVector q12_cartesian = v1_cartesian - v2_cartesian;
    float beta = (v1 + v2).Beta();
    float beta_x = beta * std::cos((v1 + v2).Phi()) * std::sin((v1 + v2).Theta());
    float beta_y = beta * std::sin((v1 + v2).Phi()) * std::sin((v1 + v2).Theta());
    float beta_z = beta * std::cos((v1 + v2).Theta());

    // longitudinally co-moving system (LCMS)
    ROOT::Math::Boost bst_z(0, 0, -beta_z); // Boost supports only PxPyPzEVector
    ROOT::Math::PxPyPzEVector q12_lcms = bst_z(q12_cartesian);
    ROOT::Math::XYZVector q_3d_lcms = q12_lcms.Vect(); // 3D q vector in LCMS
    float qout_lcms = q_3d_lcms.Dot(uv_out);
    float qside_lcms = q_3d_lcms.Dot(uv_side);
    float qlong_lcms = q_3d_lcms.Dot(uv_long);
    float qabs_lcms = q_3d_lcms.R();

    // float qabs_lcms_tmp = std::sqrt(std::pow(qout_lcms, 2) + std::pow(qside_lcms, 2) + std::pow(qlong_lcms, 2));
    // LOGF(info, "qabs_lcms = %f, qabs_lcms_tmp = %f", qabs_lcms, qabs_lcms_tmp);

    // pair rest frame (PRF)
    ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-beta_x, -beta_y, -beta_z);
    ROOT::Math::PxPyPzEVector v1_prf = boostPRF(v1_cartesian);
    ROOT::Math::PxPyPzEVector v2_prf = boostPRF(v2_cartesian);
    ROOT::Math::PxPyPzEVector rel_k = v1_prf - v2_prf;
    float kstar = 0.5 * rel_k.P();
    // LOGF(info, "qabs_lcms = %f, qinv = %f, kstar = %f", qabs_lcms, qinv, kstar);

    // ROOT::Math::PxPyPzEVector v1_lcms_cartesian = bst_z(v1_cartesian);
    // ROOT::Math::PxPyPzEVector v2_lcms_cartesian = bst_z(v2_cartesian);
    // ROOT::Math::PxPyPzEVector q12_lcms_cartesian = bst_z(q12_cartesian);
    // LOGF(info, "q12.Pz() = %f, q12_cartesian.Pz() = %f", q12.Pz(), q12_cartesian.Pz());
    // LOGF(info, "v1.Pz() = %f, v2.Pz() = %f", v1.Pz(), v2.Pz());
    // LOGF(info, "v1_lcms_cartesian.Pz() = %f, v2_lcms_cartesian.Pz() = %f", v1_lcms_cartesian.Pz(), v2_lcms_cartesian.Pz());
    // LOGF(info, "q12_lcms_cartesian.Pz() = %f", q12_lcms_cartesian.Pz());
    // LOGF(info, "q_3d_lcms.Dot(uv_out) = %f, q_3d_lcms.Dot(uv_side) = %f, q_3d.Dot(uv_out) = %f, q_3d.Dot(uv_side) = %f", q_3d_lcms.Dot(uv_out), q_3d_lcms.Dot(uv_side), q_3d.Dot(uv_out), q_3d.Dot(uv_side));
    // LOGF(info, "q12_lcms.Pz() = %f, q_3d_lcms.Dot(uv_long) = %f", q12_lcms.Pz(), q_3d_lcms.Dot(uv_long));
    // ROOT::Math::PxPyPzEVector q12_lcms_tmp = bst_z(v1_cartesian) - bst_z(v2_cartesian);
    // LOGF(info, "q12_lcms.Px() = %f, q12_lcms.Py() = %f, q12_lcms.Pz() = %f, q12_lcms_tmp.Px() = %f, q12_lcms_tmp.Py() = %f, q12_lcms_tmp.Pz() = %f", q12_lcms.Px(), q12_lcms.Py(), q12_lcms.Pz(), q12_lcms_tmp.Px(), q12_lcms_tmp.Py(), q12_lcms_tmp.Pz());
    // float qabs_lcms_tmp = q12_lcms.P();
    // LOGF(info, "qabs_lcms = %f, qabs_lcms_tmp = %f", qabs_lcms, qabs_lcms_tmp);

    if (cfgDo3D) {
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_3d"), fabs(qout_lcms), fabs(qside_lcms), fabs(qlong_lcms), kt, weight); // qosl can be [-inf, +inf] and CF is symmetric for pos and neg qosl. To reduce stat. unc. absolute value is taken here.
    } else {
      if constexpr (pairtype == ggHBTPairType::kPCMPCM) { // identical particle femtoscopy
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_1d"), qinv, qabs_lcms, kt, weight);
      } else {
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_1d"), kstar, qabs_lcms, kt, weight);
      }
    }
  }

  template <bool isTriggerAnalysis, typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2>
  void runPairing(TCollisions const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TSubInfos1 const&, TSubInfos2 const&, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCut1 const& cut1, TCut2 const& cut2)
  {
    for (auto& collision : collisions) {
      initCCDB<isTriggerAnalysis>(collision);
      int ndiphoton = 0;
      const float centralities[4] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      if constexpr (isTriggerAnalysis) {
        if (!collision.swtalias_bit(o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value))) {
          continue;
        }
        // if (collision.spherocity_ptunweighted() < cfgSpherocityMin || cfgSpherocityMax < collision.spherocity_ptunweighted()) {
        //   continue;
        // }
        // fRegistry.fill(HIST("Event/after/hSpherocity"), collision.spherocity_ptunweighted());
      }
      const float eventplanes_2_for_mix[6] = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      float ep2 = eventplanes_2_for_mix[cfgEP2Estimator_for_Mix];
      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted
      fRegistry.fill(HIST("Event/after/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

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
      std::pair<int, int64_t> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
        auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_coll = photons2.sliceBy(perCollision2, collision.globalIndex());
        for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }

          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          auto pos2 = g2.template posTrack_as<TSubInfos2>();
          auto ele2 = g2.template negTrack_as<TSubInfos2>();
          if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) { // never happens. only for protection.
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);

          float dz = g1.vz() - g2.vz();
          float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2));
          if (ggpaircuts.applydR && std::pow(dz / ggpaircuts.cfgMinDeltaZ, 2) + std::pow(dr / ggpaircuts.cfgMinDeltaR, 2) < 1.f) {
            continue;
          }
          fRegistry.fill(HIST("Pair/same/hDeltaRDeltaZ"), dr, dz, 1.f);

          fillPairHistogram<0>(collision, v1, v2, 1.f);
          ndiphoton++;

          std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
          std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, g2.globalIndex());
          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
            EMTrack g1tmp = EMTrack(ndf, g1.globalIndex(), collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0);
            g1tmp.setConversionPointXYZ(g1.vx(), g1.vy(), g1.vz());
            emh1->AddTrackToEventPool(key_df_collision, g1tmp);
            used_photonIds.emplace_back(pair_tmp_id1);
          }
          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id2) == used_photonIds.end()) {
            EMTrack g2tmp = EMTrack(ndf, g2.globalIndex(), collision.globalIndex(), g2.globalIndex(), g2.pt(), g2.eta(), g2.phi(), 0);
            g2tmp.setConversionPointXYZ(g2.vx(), g2.vy(), g2.vz());
            emh1->AddTrackToEventPool(key_df_collision, g2tmp);
            used_photonIds.emplace_back(pair_tmp_id2);
          }
        } // end of pairing loop
      } else if constexpr (pairtype == ggHBTPairType::kEEEE) {
        auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> used_pairs_per_collision;
        used_pairs_per_collision.reserve(std::pow(positrons_per_collision.size() * electrons_per_collision.size(), 2));

        for (auto& [pos1, ele1] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
          if (pos1.trackId() == ele1.trackId()) { // this is protection against pairing identical 2 tracks. // never happens. only for protection.
            continue;
          }
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut1.template IsSelectedTrack<true>(pos1, collision) || !cut1.template IsSelectedTrack<true>(ele1, collision)) {
              continue;
            }
          } else { // cut-based
            if (!cut1.template IsSelectedTrack<false>(pos1, collision) || !cut1.template IsSelectedTrack<false>(ele1, collision)) {
              continue;
            }
          }
          if (!cut1.IsSelectedPair(pos1, ele1, d_bz)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v_pos1(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v_ele1(ele1.pt(), ele1.eta(), ele1.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v1_ee = v_pos1 + v_ele1;
          float dca_pos1_3d = dca3DinSigma(pos1);
          float dca_ele1_3d = dca3DinSigma(ele1);
          float dca1_3d = std::sqrt((dca_pos1_3d * dca_pos1_3d + dca_ele1_3d * dca_ele1_3d) / 2.);
          float weight1 = 1.f;
          if (cfgApplyWeightTTCA) {
            weight1 = map_weight[std::make_pair(pos1.globalIndex(), ele1.globalIndex())];
          }

          for (auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
            if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks. // never happens. only for protection.
              continue;
            }
            if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
              if (!cut2.template IsSelectedTrack<true>(pos2, collision) || !cut2.template IsSelectedTrack<true>(ele2, collision)) {
                continue;
              }
            } else { // cut-based
              if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
                continue;
              }
            }
            if (!cut2.IsSelectedPair(pos2, ele2, d_bz)) {
              continue;
            }

            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) { // this comparison is valid in the same collision.
              continue;
            }

            float weight2 = 1.f;
            if (cfgApplyWeightTTCA) {
              weight2 = map_weight[std::make_pair(pos2.globalIndex(), ele2.globalIndex())];
            }

            float dca_pos2_3d = dca3DinSigma(pos2);
            float dca_ele2_3d = dca3DinSigma(ele2);
            float dca2_3d = std::sqrt((dca_pos2_3d * dca_pos2_3d + dca_ele2_3d * dca_ele2_3d) / 2.);

            ROOT::Math::PtEtaPhiMVector v_pos2(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele2(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v2_ee = v_pos2 + v_ele2;

            std::pair pair_tmp = std::make_pair(std::make_pair(pos1.trackId(), ele1.trackId()), std::make_pair(pos2.trackId(), ele2.trackId()));
            if (std::find(used_pairs_per_collision.begin(), used_pairs_per_collision.end(), pair_tmp) == used_pairs_per_collision.end()) {
              float deta_pos = v_pos1.Eta() - v_pos2.Eta();
              float dphi_pos = v_pos1.Phi() - v_pos2.Phi();
              o2::math_utils::bringToPMPi(dphi_pos);
              float deta_ele = v_ele1.Eta() - v_ele2.Eta();
              float dphi_ele = v_ele1.Phi() - v_ele2.Phi();
              o2::math_utils::bringToPMPi(dphi_ele);
              if (ggpaircuts.applydR && std::pow(deta_pos / ggpaircuts.cfgMinDeltaEta, 2) + std::pow(dphi_pos / ggpaircuts.cfgMinDeltaPhi, 2) < 1.f) {
                continue;
              }
              if (ggpaircuts.applydR && std::pow(deta_ele / ggpaircuts.cfgMinDeltaEta, 2) + std::pow(dphi_ele / ggpaircuts.cfgMinDeltaPhi, 2) < 1.f) {
                continue;
              }

              fRegistry.fill(HIST("Pair/same/hDeltaEtaDeltaPhi"), v_pos1.Phi() - v_pos2.Phi(), v_pos1.Eta() - v_pos2.Eta(), weight1 * weight2); // distance between 2 LS tracks
              fRegistry.fill(HIST("Pair/same/hDeltaEtaDeltaPhi"), v_ele1.Phi() - v_ele2.Phi(), v_ele1.Eta() - v_ele2.Eta(), weight1 * weight2); // distance between 2 LS tracks
              fillPairHistogram<0>(collision, v1_ee, v2_ee, weight1 * weight2);
              ndiphoton++;
              used_pairs_per_collision.emplace_back(std::make_pair(pair_tmp.first, pair_tmp.second));
              used_pairs_per_collision.emplace_back(std::make_pair(pair_tmp.second, pair_tmp.first));

              std::tuple<int, int, int, int> tuple_tmp_id1 = std::make_tuple(ndf, collision.globalIndex(), pos1.globalIndex(), ele1.globalIndex());
              std::tuple<int, int, int, int> tuple_tmp_id2 = std::make_tuple(ndf, collision.globalIndex(), pos2.globalIndex(), ele2.globalIndex());
              if (std::find(used_dileptonIds.begin(), used_dileptonIds.end(), tuple_tmp_id1) == used_dileptonIds.end()) {
                std::vector<int> possibleIds_pos1;
                std::vector<int> possibleIds_ele1;
                std::copy(pos1.ambiguousElectronsIds().begin(), pos1.ambiguousElectronsIds().end(), std::back_inserter(possibleIds_pos1));
                std::copy(ele1.ambiguousElectronsIds().begin(), ele1.ambiguousElectronsIds().end(), std::back_inserter(possibleIds_ele1));

                EMTrack g1pair = EMTrack(ndf, -1, collision.globalIndex(), -1, v1_ee.Pt(), v1_ee.Eta(), v1_ee.Phi(), v1_ee.M());
                g1pair.setGlobalPosId(pos1.globalIndex());
                g1pair.setGlobalNegId(ele1.globalIndex());
                g1pair.setPairDca3DinSigmaOTF(dca1_3d);
                g1pair.setPositiveLegPtEtaPhiM(v_pos1.Pt(), v_pos1.Eta(), v_pos1.Phi(), o2::constants::physics::MassElectron);
                g1pair.setNegativeLegPtEtaPhiM(v_ele1.Pt(), v_ele1.Eta(), v_ele1.Phi(), o2::constants::physics::MassElectron);
                g1pair.setAmbPosLegSelfIds(possibleIds_pos1);
                g1pair.setAmbNegLegSelfIds(possibleIds_ele1);
                emh1->AddTrackToEventPool(key_df_collision, g1pair);
                used_dileptonIds.emplace_back(tuple_tmp_id1);
              }
              if (std::find(used_dileptonIds.begin(), used_dileptonIds.end(), tuple_tmp_id2) == used_dileptonIds.end()) {
                std::vector<int> possibleIds_pos2;
                std::vector<int> possibleIds_ele2;
                std::copy(pos2.ambiguousElectronsIds().begin(), pos2.ambiguousElectronsIds().end(), std::back_inserter(possibleIds_pos2));
                std::copy(ele2.ambiguousElectronsIds().begin(), ele2.ambiguousElectronsIds().end(), std::back_inserter(possibleIds_ele2));

                EMTrack g2pair = EMTrack(ndf, -1, collision.globalIndex(), -1, v2_ee.Pt(), v2_ee.Eta(), v2_ee.Phi(), v2_ee.M());
                g2pair.setGlobalPosId(pos2.globalIndex());
                g2pair.setGlobalNegId(ele2.globalIndex());
                g2pair.setPairDca3DinSigmaOTF(dca2_3d);
                g2pair.setPositiveLegPtEtaPhiM(v_pos2.Pt(), v_pos2.Eta(), v_pos2.Phi(), o2::constants::physics::MassElectron);
                g2pair.setNegativeLegPtEtaPhiM(v_ele2.Pt(), v_ele2.Eta(), v_ele2.Phi(), o2::constants::physics::MassElectron);
                g2pair.setAmbPosLegSelfIds(possibleIds_pos2);
                g2pair.setAmbNegLegSelfIds(possibleIds_ele2);
                emh1->AddTrackToEventPool(key_df_collision, g2pair);
                used_dileptonIds.emplace_back(tuple_tmp_id2);
              }
            }
          } // end of g2 loop
        } // end of g1 loop
        used_pairs_per_collision.clear();
        used_pairs_per_collision.shrink_to_fit();
      } else if constexpr (pairtype == ggHBTPairType::kPCMEE) {
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

        for (auto& g1 : photons1_per_collision) {
          if (!cut1.template IsSelected<TSubInfos1>(g1)) {
            continue;
          }
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          ROOT::Math::PtEtaPhiMVector v1_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

          for (auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
            if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks. // never happens. only for protection.
              continue;
            }
            if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
              if (!cut2.template IsSelectedTrack<true>(pos2, collision) || !cut2.template IsSelectedTrack<true>(ele2, collision)) {
                continue;
              }
            } else { // cut-based
              if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
                continue;
              }
            }
            if (!cut2.IsSelectedPair(pos2, ele2, d_bz)) {
              continue;
            }

            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) { // this comparison is valid in the same collision.
              continue;
            }

            float weight = 1.f;
            if (cfgApplyWeightTTCA) {
              weight = map_weight[std::make_pair(pos2.globalIndex(), ele2.globalIndex())];
            }
            // LOGF(info, "g1.globalIndex() = %d, map_weight[std::make_pair(%d, %d)] = %f", g1.globalIndex(), pos2.globalIndex(), ele2.globalIndex(), weight);

            float dca_pos2_3d = dca3DinSigma(pos2);
            float dca_ele2_3d = dca3DinSigma(ele2);
            float dca2_3d = std::sqrt((dca_pos2_3d * dca_pos2_3d + dca_ele2_3d * dca_ele2_3d) / 2.);

            ROOT::Math::PtEtaPhiMVector v_pos2(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele2(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v2_ee = v_pos2 + v_ele2;
            float deta = v1_gamma.Eta() - v2_ee.Eta();
            float dphi = v1_gamma.Phi() - v2_ee.Phi();
            o2::math_utils::bringToPMPi(dphi);
            if (ggpaircuts.applydR && std::pow(deta / ggpaircuts.cfgMinDeltaEta, 2) + std::pow(dphi / ggpaircuts.cfgMinDeltaPhi, 2) < 1.f) {
              continue;
            }
            fRegistry.fill(HIST("Pair/same/hDeltaEtaDeltaPhi"), dphi, deta, weight);

            fillPairHistogram<0>(collision, v1_gamma, v2_ee, weight);
            ndiphoton++;
            std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
            std::tuple<int, int, int, int> tuple_tmp_id2 = std::make_tuple(ndf, collision.globalIndex(), pos2.globalIndex(), ele2.globalIndex());
            if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
              emh1->AddTrackToEventPool(key_df_collision, EMTrack(ndf, g1.globalIndex(), collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0));
              used_photonIds.emplace_back(pair_tmp_id1);
            }
            if (std::find(used_dileptonIds.begin(), used_dileptonIds.end(), tuple_tmp_id2) == used_dileptonIds.end()) {
              EMTrack g2pair = EMTrack(ndf, -1, collision.globalIndex(), -1, v2_ee.Pt(), v2_ee.Eta(), v2_ee.Phi(), v2_ee.M());
              g2pair.setPairDca3DinSigmaOTF(dca2_3d);
              emh2->AddTrackToEventPool(key_df_collision, g2pair);
              used_dileptonIds.emplace_back(tuple_tmp_id2);
            }
          } // end of g2 loop
        } // end of g1 loop
      }

      // event mixing
      if (!cfgDoMix || !(ndiphoton > 0)) {
        continue;
      }

      // make a vector of selected photons in this collision.
      auto selected_photons1_in_this_event = emh1->GetTracksPerCollision(key_df_collision);
      auto selected_photons2_in_this_event = emh2->GetTracksPerCollision(key_df_collision);

      auto collisionIds1_in_mixing_pool = emh1->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIds2_in_mixing_pool = emh2->GetCollisionIdsFromEventPool(key_bin);

      if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
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

          for (auto& g1 : selected_photons1_in_this_event) {
            for (auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);

              float dz = g1.vz() - g2.vz();
              float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2));
              if (ggpaircuts.applydR && std::pow(dz / ggpaircuts.cfgMinDeltaZ, 2) + std::pow(dr / ggpaircuts.cfgMinDeltaR, 2) < 1.f) {
                continue;
              }
              fRegistry.fill(HIST("Pair/mix/hDeltaRDeltaZ"), dr, dz, 1.f);
              fillPairHistogram<1>(collision, v1, v2, 1.f);
            }
          }
        } // end of loop over mixed event pool
      } else if constexpr (pairtype == ggHBTPairType::kEEEE) {
        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
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

          for (auto& g1 : selected_photons1_in_this_event) {
            for (auto& g2 : photons1_from_event_pool) {
              auto v1amb_pos_Ids = g1.ambiguousPosLegIds();
              auto v1amb_neg_Ids = g1.ambiguousNegLegIds();
              auto v2amb_pos_Ids = g2.ambiguousPosLegIds();
              auto v2amb_neg_Ids = g2.ambiguousNegLegIds();

              bool is_found_pos1 = std::find(v2amb_pos_Ids.begin(), v2amb_pos_Ids.end(), g1.globalIndexPos()) != v2amb_pos_Ids.end();
              bool is_found_neg1 = std::find(v2amb_neg_Ids.begin(), v2amb_neg_Ids.end(), g1.globalIndexPos()) != v2amb_neg_Ids.end();
              bool is_found_pos2 = std::find(v1amb_pos_Ids.begin(), v1amb_pos_Ids.end(), g2.globalIndexPos()) != v1amb_pos_Ids.end();
              bool is_found_neg2 = std::find(v1amb_neg_Ids.begin(), v1amb_neg_Ids.end(), g2.globalIndexPos()) != v1amb_neg_Ids.end();
              // LOGF(info, "is_found_pos1 = %d, is_found_neg1 = %d, is_found_pos2 = %d, is_found_neg2 = %d", is_found_pos1, is_found_neg1, is_found_pos2, is_found_neg2);

              auto pos1 = g1.getPositiveLeg();
              auto ele1 = g1.getNegativeLeg();
              auto pos2 = g2.getPositiveLeg();
              auto ele2 = g2.getNegativeLeg();

              if ((g1.dfId() == g2.dfId()) && ((is_found_pos1 && is_found_pos2) || (is_found_neg1 && is_found_neg2))) {
                // LOGF(info, "event id = %d: same track is found. t1.globalIndex() = %d, t1.sign() = %d, t1.pt() = %f, t1.eta() = %f, t1.phi() = %f, t2.globalIndex() = %d, t2.sign() = %d, t2.pt() = %f, t2.eta() = %f, t2.phi() = %f, deta = %f, dphi = %f (rad.)", ev_id, t1.globalIndex(), t1.sign(), t1.pt(), t1.eta(), t1.phi(), t2.globalIndex(), t2.sign(), t2.pt(), t2.eta(), t2.phi(), t1.eta() - t2.eta(), t1.phi() - t2.phi());
                continue; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
              }

              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), g1.mass());
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
              float deta_pos = pos1.Eta() - pos2.Eta();
              float dphi_pos = pos1.Phi() - pos2.Phi();
              o2::math_utils::bringToPMPi(dphi_pos);
              float deta_ele = ele1.Eta() - ele2.Eta();
              float dphi_ele = ele1.Phi() - ele2.Phi();
              o2::math_utils::bringToPMPi(dphi_ele);
              if (ggpaircuts.applydR && std::pow(deta_pos / ggpaircuts.cfgMinDeltaEta, 2) + std::pow(dphi_pos / ggpaircuts.cfgMinDeltaPhi, 2) < 1.f) {
                continue;
              }
              if (ggpaircuts.applydR && std::pow(deta_ele / ggpaircuts.cfgMinDeltaEta, 2) + std::pow(dphi_ele / ggpaircuts.cfgMinDeltaPhi, 2) < 1.f) {
                continue;
              }
              fRegistry.fill(HIST("Pair/mix/hDeltaEtaDeltaPhi"), pos1.Phi() - pos2.Phi(), pos1.Eta() - pos2.Eta(), 1.f);
              fRegistry.fill(HIST("Pair/mix/hDeltaEtaDeltaPhi"), ele1.Phi() - ele2.Phi(), ele1.Eta() - ele2.Eta(), 1.f);
              fillPairHistogram<1>(collision, v1, v2, 1.f);
            }
          }
        } // end of loop over mixed event pool
      } else if constexpr (pairtype == ggHBTPairType::kPCMEE) { // [photon1 from event1, photon2 from event2] and [photon1 from event2, photon2 from event1]
        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
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

          for (auto& g1 : selected_photons1_in_this_event) {                    // PCM
            for (auto& g2 : photons2_from_event_pool) {                         // dielectron
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.0); // keep v1 for PCM
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
              float deta = v1.Eta() - v2.Eta();
              float dphi = v1.Phi() - v2.Phi();
              o2::math_utils::bringToPMPi(dphi);
              if (ggpaircuts.applydR && std::pow(deta / ggpaircuts.cfgMinDeltaEta, 2) + std::pow(dphi / ggpaircuts.cfgMinDeltaPhi, 2) < 1.f) {
                continue;
              }
              fRegistry.fill(HIST("Pair/mix/hDeltaEtaDeltaPhi"), dphi, deta, 1.f);
              fillPairHistogram<1>(collision, v1, v2, 1.f);
            }
          }
        } // end of loop over mixed event pool2

        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
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

          for (auto& g1 : selected_photons2_in_this_event) {                    // dielectron
            for (auto& g2 : photons1_from_event_pool) {                         // PCM
              ROOT::Math::PtEtaPhiMVector v1(g2.pt(), g2.eta(), g2.phi(), 0.0); // keep v1 for PCM
              ROOT::Math::PtEtaPhiMVector v2(g1.pt(), g1.eta(), g1.phi(), g1.mass());
              float deta = v1.Eta() - v2.Eta();
              float dphi = v1.Phi() - v2.Phi();
              o2::math_utils::bringToPMPi(dphi);
              if (ggpaircuts.applydR && std::pow(deta / ggpaircuts.cfgMinDeltaEta, 2) + std::pow(dphi / ggpaircuts.cfgMinDeltaPhi, 2) < 1.f) {
                continue;
              }
              fRegistry.fill(HIST("Pair/mix/hDeltaEtaDeltaPhi"), dphi, deta, 1.f);
              fillPairHistogram<1>(collision, v1, v2, 1.f);
            }
          }
        } // end of loop over mixed event pool1
      }

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh2->AddCollisionIdAtLast(key_bin, key_df_collision);
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
      }
    } // end of collision loop
  }

  std::map<std::pair<int, int>, float> map_weight; // <posId, negId> -> float
  template <bool isTriggerAnalysis, typename TCollisions, typename TTracks, typename TCut>
  void fillDileptonPairWeightMap(TCollisions const& collisions, TTracks const& tracks, TCut const& cut)
  {
    std::vector<std::pair<int, int>> passed_pairIds;
    passed_pairIds.reserve(positrons.size() * electrons.size());

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
        // if (collision.spherocity_ptunweighted() < cfgSpherocityMin || cfgSpherocityMax < collision.spherocity_ptunweighted()) {
        //   continue;
        // }
      }
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
        if (pos.trackId() == ele.trackId()) { // this is protection against pairing identical 2 tracks. // never happens. only for protection.
          continue;
        }
        if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
          if (!cut.template IsSelectedTrack<true>(pos, collision) || !cut.template IsSelectedTrack<true>(ele, collision)) {
            continue;
          }
        } else { // cut-based
          if (!cut.template IsSelectedTrack<false>(pos, collision) || !cut.template IsSelectedTrack<false>(ele, collision)) {
            continue;
          }
        }
        if (!cut.IsSelectedPair(pos, ele, d_bz)) {
          continue;
        }
        passed_pairIds.emplace_back(std::make_pair(pos.globalIndex(), ele.globalIndex()));
      } // end of dielectron pairing loop
    } // end of collision loop

    for (auto& pairId : passed_pairIds) {
      auto t1 = tracks.rawIteratorAt(std::get<0>(pairId));
      auto t2 = tracks.rawIteratorAt(std::get<1>(pairId));
      // LOGF(info, "std::get<0>(pairId) = %d, std::get<1>(pairId) = %d, t1.globalIndex() = %d, t2.globalIndex() = %d", std::get<0>(pairId), std::get<1>(pairId), t1.globalIndex(), t2.globalIndex());

      float n = 1.f; // include myself.
      for (auto& ambId1 : t1.ambiguousElectronsIds()) {
        for (auto& ambId2 : t2.ambiguousElectronsIds()) {
          if (std::find(passed_pairIds.begin(), passed_pairIds.end(), std::make_pair(ambId1, ambId2)) != passed_pairIds.end()) {
            // LOGF(info, "repeated pair is found. t1.globalIndex() = %d, t2.globalIndex() = %d, ambId1 = %d, ambId2 = %d", t1.globalIndex(), t2.globalIndex(), ambId1, ambId2);
            n += 1.f;
          }
        }
      }
      map_weight[pairId] = 1.f / n;

    } // end of passed_pairIds loop

    passed_pairIds.clear();
    passed_pairIds.shrink_to_fit();
  }

  Filter trackFilter = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter = (dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < dielectroncuts.cfg_min_TPCNsigmaPi || dielectroncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi);
  Filter ttcaFilter = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);

  Partition<FilteredMyTracks> positrons = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyTracks> electrons = o2::aod::emprimaryelectron::sign < int8_t(0);

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrack>;
  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;
  std::vector<std::pair<int, int>> used_photonIds;              // <ndf, trackId>
  std::vector<std::tuple<int, int, int, int>> used_dileptonIds; // <ndf, trackId>
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;
  Preslice<FilteredMyTracks> perCollision_electron = aod::emprimaryelectron::emeventId;

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_multiplicity = cfgNtracksPV08Min <= o2::aod::mult::multNTracksPV && o2::aod::mult::multNTracksPV < cfgNtracksPV08Max;
  Filter collisionFilter_occupancy = eventcuts.cfgOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgOccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  int ndf = 0;
  void processAnalysis(FilteredMyCollisions const& collisions, Types const&... args)
  {
    if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      runPairing<false>(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut);
    } else if constexpr (pairtype == ggHBTPairType::kPCMEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillDileptonPairWeightMap<false>(collisions, emprimaryelectrons, fDielectronCut);
      }
      runPairing<false>(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDielectronCut);
    } else if constexpr (pairtype == ggHBTPairType::kEEEE) {
      auto emprimaryelectrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillDileptonPairWeightMap<false>(collisions, emprimaryelectrons, fDielectronCut);
      }
      runPairing<false>(collisions, nullptr, nullptr, emprimaryelectrons, emprimaryelectrons, perCollision_electron, perCollision_electron, fDielectronCut, fDielectronCut);
    }
    ndf++;
  }
  PROCESS_SWITCH(PhotonHBT, processAnalysis, "pairing for analysis", false);

  using FilteredMyCollisionsWithSWT = soa::Filtered<MyCollisionsWithSWT>;
  void processTriggerAnalysis(FilteredMyCollisionsWithSWT const& collisions, Types const&... args)
  {
    if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      runPairing<true>(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut);
    } else if constexpr (pairtype == ggHBTPairType::kPCMEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillDileptonPairWeightMap<true>(collisions, emprimaryelectrons, fDielectronCut);
      }
      runPairing<true>(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDielectronCut);
    } else if constexpr (pairtype == ggHBTPairType::kEEEE) {
      auto emprimaryelectrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillDileptonPairWeightMap<true>(collisions, emprimaryelectrons, fDielectronCut);
      }
      runPairing<true>(collisions, nullptr, nullptr, emprimaryelectrons, emprimaryelectrons, perCollision_electron, perCollision_electron, fDielectronCut, fDielectronCut);
    }
    ndf++;
  }
  PROCESS_SWITCH(PhotonHBT, processTriggerAnalysis, "pairing analysis on trigger data", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(PhotonHBT, processDummy, "Dummy function", true);
};

#endif // PWGEM_DILEPTON_CORE_PHOTONHBT_H_
