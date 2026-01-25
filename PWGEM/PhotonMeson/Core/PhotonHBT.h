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

/// \file PhotonHBT.h
/// \brief This code loops over v0 photons and makes pairs for photon HBT analysis.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_PHOTONHBT_H_
#define PWGEM_PHOTONMESON_CORE_PHOTONHBT_H_

#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
//
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h> // IWYU pragma: keep
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace o2::aod::pwgem::photon::core::photonhbt
{
enum class ggHBTPairType : int {
  kPCMPCM = 0,
};
} // namespace o2::aod::pwgem::photon::core::photonhbt

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils;
using namespace o2::aod::pwgem::photon::core::photonhbt;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

template <ggHBTPairType pairtype, typename... Types>
struct PhotonHBT {
  // Configurables

  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  // Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<bool> cfgDo3D{"cfgDo3D", false, "enable 3D analysis"};
  Configurable<int> cfgEP2Estimator_for_Mix{"cfgEP2Estimator_for_Mix", 3, "FT0M:0, FT0A:1, FT0C:2, BTot:3, BPos:4, BNeg:5"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.8, "maximum rapidity for reconstructed particles"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {16, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
  Configurable<bool> cfgUseLCMS{"cfgUseLCMS", true, "measure relative momentum in LCMS for 1D"}; // always in LCMS for 3D

  ConfigurableAxis ConfQBins{"ConfQBins", {60, 0, +0.3f}, "q bins for output histograms"};
  ConfigurableAxis ConfKtBins{"ConfKtBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}, "kT bins for output histograms"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
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
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 16.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.997, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
    Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};

    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
  } pcmcuts;

  struct : ConfigurableGroup {
    std::string prefix = "ggpaircut_group";
    Configurable<float> cfgMinDR_CosOA{"cfgMinDR_CosOA", -1, "min. dr/cosOA for kPCMPCM"};
  } ggpaircuts;

  ~PhotonHBT()
  {
    delete emh1;
    emh1 = 0x0;
    delete emh2;
    emh2 = 0x0;

    map_mixed_eventId_to_globalBC.clear();

    used_photonIds_per_col.clear();
    used_photonIds_per_col.shrink_to_fit();
  }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  // static constexpr std::string_view event_types[2] = {"before", "after"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::mt19937 engine;
  std::uniform_int_distribution<int> dist01;

  // o2::ccdb::CcdbApi ccdbApi;
  // Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  // float d_bz;

  std::vector<float> zvtx_bin_edges;
  std::vector<float> cent_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;

  void init(InitContext& /*context*/)
  {
    mRunNumber = 0;
    // d_bz = 0;

    // ccdb->setURL(ccdburl);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setFatalWhenNull(false);

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

    emh1 = new MyEMH(ndepth);
    emh2 = new MyEMH(ndepth);

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    DefineEMEventCut();
    DefinePCMCut();

    addhistograms();

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_int_distribution<int>(0, 1);

    fRegistry.add("Pair/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{10001, -0.5, 10000.5}}, true);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    // // In case override, don't proceed, please - no CCDB access required
    // if (d_bz_input > -990) {
    //   d_bz = d_bz_input;
    //   o2::parameters::GRPMagField grpmag;
    //   if (std::fabs(d_bz) > 1e-5) {
    //     grpmag.setL3Current(30000.f / (d_bz / 5.0f));
    //   }
    //   mRunNumber = collision.runNumber();
    //   return;
    // }

    // auto run3grp_timestamp = collision.timestamp();
    // o2::parameters::GRPObject* grpo = 0x0;
    // o2::parameters::GRPMagField* grpmag = 0x0;
    // if (!skipGRPOquery)
    //   grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    // if (grpo) {
    //   // Fetch magnetic field from ccdb for current collision
    //   d_bz = grpo->getNominalL3Field();
    //   LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    // } else {
    //   grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
    //   if (!grpmag) {
    //     LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
    //   }
    //   // Fetch magnetic field from ccdb for current collision
    //   d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
    //   LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    // }
    mRunNumber = collision.runNumber();
  }

  void addhistograms()
  {
    // o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);
    static constexpr std::string_view qvec_det_names[6] = {"FT0M", "FT0A", "FT0C", "BTot", "BPos", "BNeg"};
    fRegistry.add("Event/before/hEP2_CentFT0C_forMix", Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEP2Estimator_for_Mix].data()), kTH2F, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);
    fRegistry.add("Event/after/hEP2_CentFT0C_forMix", Form("2nd harmonics event plane for mix;centrality FT0C (%%);#Psi_{2}^{%s} (rad.)", qvec_det_names[cfgEP2Estimator_for_Mix].data()), kTH2F, {{110, 0, 110}, {180, -o2::constants::math::PIHalf, +o2::constants::math::PIHalf}}, false);

    // pair info
    const AxisSpec axis_kt{ConfKtBins, "k_{T} (GeV/c)"};
    const AxisSpec axis_qinv{ConfQBins, "q_{inv} (GeV/c)"};
    const AxisSpec axis_qabs_lcms{ConfQBins, "|#bf{q}|^{LCMS} (GeV/c)"};
    const AxisSpec axis_qout{ConfQBins, "q_{out} (GeV/c)"};   // qout does not change between LAB and LCMS frame
    const AxisSpec axis_qside{ConfQBins, "q_{side} (GeV/c)"}; // qside does not change between LAB and LCMS frame
    const AxisSpec axis_qlong{ConfQBins, "q_{long} (GeV/c)"};

    if (cfgDo3D) { // 3D
      fRegistry.add("Pair/same/hs_3d", "diphoton correlation 3D LCMS", kTHnSparseD, {axis_qout, axis_qside, axis_qlong, axis_kt}, true);
    } else { // 1D
      if (cfgUseLCMS) {
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D LCMS", kTHnSparseD, {axis_qabs_lcms, axis_kt}, true);
      } else {
        fRegistry.add("Pair/same/hs_1d", "diphoton correlation 1D", kTHnSparseD, {axis_qinv, axis_kt}, true);
      }
    }

    if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
      fRegistry.add("Pair/same/hDeltaRCosOA", "distance between 2 conversion points;#Deltar/cos(#theta_{op}/2) (cm)", kTH1D, {{100, 0, 100}}, true); // dr/cosOA of conversion points
    }

    fRegistry.addClone("Pair/same/", "Pair/mix/");
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
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
    fEMEventCut.SetRequireGoodITSLayer3(eventcuts.cfgRequireGoodITSLayer3);
    fEMEventCut.SetRequireGoodITSLayer0123(eventcuts.cfgRequireGoodITSLayer0123);
    fEMEventCut.SetRequireGoodITSLayersAll(eventcuts.cfgRequireGoodITSLayersAll);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, 1e10f);
    fV0PhotonCut.SetV0EtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
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
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfg_disable_itsonly_track);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfg_disable_tpconly_track);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetRequireITSTPC(pcmcuts.cfg_require_v0_with_itstpc);
    fV0PhotonCut.SetRequireITSonly(pcmcuts.cfg_require_v0_with_itsonly);
    fV0PhotonCut.SetRequireTPConly(pcmcuts.cfg_require_v0_with_tpconly);
  }

  template <int ev_id, typename TCollision>
  void fillPairHistogram(TCollision const&, const ROOT::Math::PtEtaPhiMVector v1, const ROOT::Math::PtEtaPhiMVector v2, const float weight = 1.f)
  {
    float rndm = std::pow(-1, dist01(engine) % 2); // +1 or -1 to randomize order between 1 and 2.
    // Lab. frame
    ROOT::Math::PtEtaPhiMVector q12 = (v1 - v2) * rndm;
    ROOT::Math::PtEtaPhiMVector k12 = 0.5 * (v1 + v2);
    float qinv = -q12.M(); // for identical particles -> qinv = 2 x kstar
    float kt = k12.Pt();

    ROOT::Math::XYZVector uv_out(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0); // unit vector for out. i.e. parallel to kt
    ROOT::Math::XYZVector uv_long(0, 0, 1);                                    // unit vector for long, beam axis
    ROOT::Math::XYZVector uv_side = uv_out.Cross(uv_long);                     // unit vector for side

    ROOT::Math::PxPyPzEVector v1_cartesian(v1);
    ROOT::Math::PxPyPzEVector v2_cartesian(v2);
    ROOT::Math::PxPyPzEVector q12_cartesian = (v1_cartesian - v2_cartesian) * rndm;
    float beta = (v1 + v2).Beta();
    // float beta_x = beta * std::cos((v1 + v2).Phi()) * std::sin((v1 + v2).Theta());
    // float beta_y = beta * std::sin((v1 + v2).Phi()) * std::sin((v1 + v2).Theta());
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

    // // pair rest frame (PRF)
    // ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-beta_x, -beta_y, -beta_z);
    // ROOT::Math::PxPyPzEVector v1_prf = boostPRF(v1_cartesian);
    // ROOT::Math::PxPyPzEVector v2_prf = boostPRF(v2_cartesian);
    // ROOT::Math::PxPyPzEVector rel_k = (v1_prf - v2_prf) * rndm;
    // float kstar = 0.5 * rel_k.P();
    // // LOGF(info, "qabs_lcms = %f, qinv = %f, kstar = %f", qabs_lcms, qinv, kstar);

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
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_3d"), std::fabs(qout_lcms), std::fabs(qside_lcms), std::fabs(qlong_lcms), kt, weight); // qosl can be [-inf, +inf] and CF is symmetric for pos and neg qosl. To reduce stat. unc. absolute value is taken here.
    } else {
      if (cfgUseLCMS) {
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_1d"), qabs_lcms, kt, weight);
      } else {
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_1d"), qinv, kt, weight);
      }
    }
  }

  template <typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2>
  void runPairing(TCollisions const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TSubInfos1 const&, TSubInfos2 const&, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCut1 const& cut1, TCut2 const& cut2)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      int ndiphoton = 0;
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      const float eventplanes_2_for_mix[6] = {collision.ep2ft0m(), collision.ep2ft0a(), collision.ep2ft0c(), collision.ep2btot(), collision.ep2bpos(), collision.ep2bneg()};
      float ep2 = eventplanes_2_for_mix[cfgEP2Estimator_for_Mix];
      fRegistry.fill(HIST("Event/before/hEP2_CentFT0C_forMix"), collision.centFT0C(), ep2);

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
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

      auto key_bin = std::make_tuple(zbin, centbin, epbin, occbin);
      auto key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
        auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_coll = photons2.sliceBy(perCollision2, collision.globalIndex());
        for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {
          if (!cut1.template IsSelected<decltype(g1), TSubInfos1>(g1) || !cut2.template IsSelected<decltype(g2), TSubInfos2>(g2)) {
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

          float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2) + std::pow(g1.vz() - g2.vz(), 2));
          ROOT::Math::XYZVector cp1(g1.vx(), g1.vy(), g1.vz());
          ROOT::Math::XYZVector cp2(g2.vx(), g2.vy(), g2.vz());
          float opa = std::acos(cp1.Dot(cp2) / (std::sqrt(cp1.Mag2()) * std::sqrt(cp2.Mag2()))); // opening angle between 2 conversion points
          o2::math_utils::bringTo02Pi(opa);
          if (opa > o2::constants::math::PI) {
            opa -= o2::constants::math::PI;
          }
          float cosOA = std::cos(opa / 2.f);
          if (dr / cosOA < ggpaircuts.cfgMinDR_CosOA) {
            continue;
          }
          fRegistry.fill(HIST("Pair/same/hDeltaRCosOA"), dr / cosOA);

          fillPairHistogram<0>(collision, v1, v2, 1.f);
          ndiphoton++;

          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g1.globalIndex()) == used_photonIds_per_col.end()) {
            EMPair g1tmp = EMPair(g1.pt(), g1.eta(), g1.phi(), 0);
            g1tmp.setConversionPointXYZ(g1.vx(), g1.vy(), g1.vz());
            emh1->AddTrackToEventPool(key_df_collision, g1tmp);
            used_photonIds_per_col.emplace_back(g1.globalIndex());
          }
          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g2.globalIndex()) == used_photonIds_per_col.end()) {
            EMPair g2tmp = EMPair(g2.pt(), g2.eta(), g2.phi(), 0);
            g2tmp.setConversionPointXYZ(g2.vx(), g2.vy(), g2.vz());
            emh1->AddTrackToEventPool(key_df_collision, g2tmp);
            used_photonIds_per_col.emplace_back(g2.globalIndex());
          }
        } // end of pairing loop
      }

      used_photonIds_per_col.clear();
      used_photonIds_per_col.shrink_to_fit();

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

              float dr = std::sqrt(std::pow(g1.vx() - g2.vx(), 2) + std::pow(g1.vy() - g2.vy(), 2) + std::pow(g1.vz() - g2.vz(), 2));
              ROOT::Math::XYZVector cp1(g1.vx(), g1.vy(), g1.vz());
              ROOT::Math::XYZVector cp2(g2.vx(), g2.vy(), g2.vz());
              float opa = std::acos(cp1.Dot(cp2) / (std::sqrt(cp1.Mag2()) * std::sqrt(cp2.Mag2()))); // opening angle between 2 conversion points
              o2::math_utils::bringTo02Pi(opa);
              if (opa > o2::constants::math::PI) {
                opa -= o2::constants::math::PI;
              }
              float cosOA = std::cos(opa / 2.f);
              if (dr / cosOA < ggpaircuts.cfgMinDR_CosOA) {
                continue;
              }
              fRegistry.fill(HIST("Pair/mix/hDeltaRCosOA"), dr / cosOA);

              fillPairHistogram<1>(collision, v1, v2, 1.f);
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

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMPair>;
  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;
  std::vector<int> used_photonIds_per_col; // <trackId>
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  int ndf = 0;
  void processAnalysis(FilteredMyCollisions const& collisions, Types const&... args)
  {
    if constexpr (pairtype == ggHBTPairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      runPairing(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut);
    }
    ndf++;
  }
  PROCESS_SWITCH(PhotonHBT, processAnalysis, "pairing for analysis", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(PhotonHBT, processDummy, "Dummy function", true);
};

#endif // PWGEM_PHOTONMESON_CORE_PHOTONHBT_H_
