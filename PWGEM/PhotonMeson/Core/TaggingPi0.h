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

/// \file TaggingPi0.h
/// \brief This code loops over photons and makes pairs for direct photon analysis.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_TAGGINGPI0_H_
#define PWGEM_PHOTONMESON_CORE_TAGGINGPI0_H_

#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::photonpair;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithJJMC = soa::Join<MyCollisions, aod::EMEventsWeight>;
using MyCollisionWithJJMC = MyCollisionsWithJJMC::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyPrimaryElectrons = soa::Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronEMEventIds>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMCEMEventIds>;
using MyEMCCluster = MyEMCClusters::iterator;

using MyPHOSClusters = soa::Join<aod::PHOSClusters, aod::PHOSEMEventIds>;
using MyPHOSCluster = MyPHOSClusters::iterator;

template <PairType pairtype, typename... Types>
struct TaggingPi0 {
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};

  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {VARIABLE_WIDTH, -o2::constants::math::PIHalf, -o2::constants::math::PIQuarter, 0.0f, +o2::constants::math::PIQuarter, +o2::constants::math::PIHalf}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};
  ConfigurableAxis ConfPtBins{"ConfPtBins", {100, 0, 10}, "pT bins for output histograms"};

  EMPhotonEventCut fEMEventCut;
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
    Configurable<bool> cfgRequireEMCReadoutInMB{"cfgRequireEMCReadoutInMB", false, "require the EMC to be read out in an MB collision (kTVXinEMC)"};
    Configurable<bool> cfgRequireEMCHardwareTriggered{"cfgRequireEMCHardwareTriggered", false, "require the EMC to be hardware triggered (kEMC7 or kDMC7)"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_pt_v0{"cfg_max_pt_v0", 1e+10, "max pT for v0 photons at PV"};
    Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.997, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
    Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};

    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 10, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};
  } pcmcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.1, "max mass"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.8, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.05, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.05, "max dca Z for single track in cm"};
    Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1.5, "max DCA 3D in sigma"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTOFif), "pid scheme [kTOFif : 0, kTPConly : 1]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -0.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +0.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
  } dileptoncuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccut_group";
    Configurable<std::string> clusterDefinition{"clusterDefinition", "kV3Default", "Clusterizer to be selected, e.g. V3Default"};
    Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};
    Configurable<float> EMC_minTime{"EMC_minTime", -20., "Minimum cluster time for EMCal time cut"};
    Configurable<float> EMC_maxTime{"EMC_maxTime", +25., "Maximum cluster time for EMCal time cut"};
    Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    Configurable<int> cfgDistanceToEdge{"cfgDistanceToEdge", 1, "Distance to edge in cells required for rotated cluster to be accepted"};
  } emccuts;

  PHOSPhotonCut fPHOSCut;
  struct : ConfigurableGroup {
    std::string prefix = "phoscut_group";
    Configurable<float> cfg_min_Ecluster{"cfg_min_Ecluster", 0.3, "Minimum cluster energy for PHOS in GeV"};
  } phoscuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> zvtx_bin_edges;
  std::vector<float> cent_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  void init(InitContext&)
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

    emh1 = new MyEMH(ndepth);
    emh2 = new MyEMH(ndepth);

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);

    addHistogrms();
    DefineEMEventCut();
    DefinePCMCut();
    DefineDileptonCut();
    DefineEMCCut();
    DefinePHOSCut();

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
    mRunNumber = collision.runNumber();
  }

  ~TaggingPi0()
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

  void addHistogrms()
  {
    TString mggTitle = "ee#gamma";
    if constexpr (pairtype == PairType::kPCMDalitzEE) {
      mggTitle = "ee#gamma";
    } else {
      mggTitle = "#gamma#gamma";
    }

    const AxisSpec axis_m{200, 0, 0.4, Form("m_{%s} (GeV/c^{2})", mggTitle.Data())};
    const AxisSpec axis_pt{ConfPtBins, "p_{T,#gamma} (GeV/c)"};

    fRegistry.add("Photon/hPt", "p_{T,#gamma};p_{T,#gamma} (GeV/c)", kTH1D, {axis_pt}, true);
    fRegistry.add("Photon/hEtaPhi", "#varphi vs. #eta;#varphi_{#gamma} (rad.);#eta_{#gamma}", kTH2D, {{90, 0, o2::constants::math::TwoPI}, {40, -1, +1}}, true);
    fRegistry.add("Pair/same/hMvsPt", "mass vs. p_{T,#gamma}", kTH2D, {axis_m, axis_pt}, true);
    fRegistry.addClone("Pair/same/", "Pair/mix/");
    fRegistry.add("Pair/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{10001, -0.5, 10000.5}}, true);
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
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
  }

  void DefinePHOSCut()
  {
    fPHOSCut.SetEnergyRange(phoscuts.cfg_min_Ecluster, 1e+10);
  }

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;
  Preslice<MyEMCClusters> perCollision_emc = aod::emccluster::emeventId;
  Preslice<MyPHOSClusters> perCollision_phos = aod::phoscluster::emeventId;

  Preslice<MyPrimaryElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Partition<MyPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);
  Partition<MyPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrack>;
  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;
  std::vector<int> used_photonIds_per_col;                   // <ndf, trackId>
  std::vector<std::pair<int, int>> used_dileptonIds_per_col; // <ndf, trackId>
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  template <typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2>
  void runPairing(TCollisions const& collisions,
                  TPhotons1 const& photons1, TPhotons2 const& photons2,
                  TSubInfos1 const&, TSubInfos2 const&,
                  TPreslice1 const& perCollision1, TPreslice2 const& perCollision2,
                  TCut1 const& cut1, TCut2 const& cut2)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      int ndiphoton = 0;

      float weight = 1.f;
      if constexpr (std::is_same_v<std::decay_t<TCollisions>, FilteredMyCollisionsWithJJMC>) {
        weight = collision.weight();
      }

      if (eventcuts.onlyKeepWeightedEvents && std::fabs(weight - 1.f) < 1e-10) {
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

      if constexpr (pairtype == PairType::kPCMDalitzEE) {
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());                                         // PCM
        auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache); // positrons
        auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache); // electrons

        for (const auto& g1 : photons1_per_collision) {
          if (!cut1.template IsSelected<decltype(g1), TSubInfos1>(g1)) {
            continue;
          }
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);
          fRegistry.fill(HIST("Photon/hPt"), v_gamma.Pt(), weight);
          fRegistry.fill(HIST("Photon/hEtaPhi"), v_gamma.Phi() > 0 ? v_gamma.Phi() : v_gamma.Phi() + o2::constants::math::TwoPI, v_gamma.Eta(), weight);

          for (const auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {

            if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks.
              continue;
            }
            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
              continue;
            }

            if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
              continue;
            }

            if (!cut2.IsSelectedPair(pos2, ele2, d_bz)) {
              continue;
            }

            ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
            ROOT::Math::PtEtaPhiMVector veeg = v_gamma + v_pos + v_ele;
            fRegistry.fill(HIST("Pair/same/hMvsPt"), veeg.M(), v_gamma.Pt(), weight);

            std::pair<int, int> tuple_tmp_id2 = std::make_pair(pos2.trackId(), ele2.trackId());
            if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g1.globalIndex()) == used_photonIds_per_col.end()) {
              emh1->AddTrackToEventPool(key_df_collision, EMTrack(g1.pt(), g1.eta(), g1.phi(), 0));
              used_photonIds_per_col.emplace_back(g1.globalIndex());
            }
            if (std::find(used_dileptonIds_per_col.begin(), used_dileptonIds_per_col.end(), tuple_tmp_id2) == used_dileptonIds_per_col.end()) {
              emh2->AddTrackToEventPool(key_df_collision, EMTrack(v_ee.Pt(), v_ee.Eta(), v_ee.Phi(), v_ee.M()));
              used_dileptonIds_per_col.emplace_back(tuple_tmp_id2);
            }
            ndiphoton++;
          } // end of dielectron loop
        } // end of g1 loop
      } else {                                                                                  // PCM-EMC, PCM-PHOS. Not supported.
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex()); // PCM
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex()); // EMC or PHOS

        for (const auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<decltype(g1), TSubInfos1>(g1) || !cut2.template IsSelected<decltype(g2), TSubInfos2>(g2)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

          fRegistry.fill(HIST("Pair/same/hMvsPt"), v12.M(), v1.Pt(), weight);

          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g1.globalIndex()) == used_photonIds_per_col.end()) {
            emh1->AddTrackToEventPool(key_df_collision, EMTrack(g1.pt(), g1.eta(), g1.phi(), 0));
            used_photonIds_per_col.emplace_back(g1.globalIndex());
          }
          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g2.globalIndex()) == used_photonIds_per_col.end()) {
            emh2->AddTrackToEventPool(key_df_collision, EMTrack(g2.pt(), g2.eta(), g2.phi(), 0));
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
      // auto selected_photons2_in_this_event = emh2->GetTracksPerCollision(key_df_collision);

      // auto collisionIds1_in_mixing_pool = emh1->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIds2_in_mixing_pool = emh2->GetCollisionIdsFromEventPool(key_bin);

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

        for (const auto& g1 : selected_photons1_in_this_event) { // [photon from event1, dilepton from event2]
          for (const auto& g2 : photons2_from_event_pool) {
            ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
            ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
            if constexpr (pairtype == PairType::kPCMDalitzEE) {
              v2.SetM(g2.mass());
            }
            ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
            fRegistry.fill(HIST("Pair/mix/hMvsPt"), v12.M(), v1.Pt(), weight);
          }
        }
      } // end of loop over mixed event pool

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh2->AddCollisionIdAtLast(key_bin, key_df_collision);
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
      }

    } // end of collision loop
  }

  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  int ndf = 0;
  void processAnalysis(FilteredMyCollisions const& collisions, Types const&... args)
  {
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == PairType::kPCMDalitzEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runPairing(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut);
    }
    // else if constexpr (pairtype == PairType::kPCMEMC) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto emcclusters = std::get<2>(std::tie(args...));
    //   auto emcmatchedtracks = std::get<3>(std::tie(args...));
    //   runPairing(collisions, v0photons, emcclusters, v0legs, nullptr, perCollision_pcm, perCollision_emc, fV0PhotonCut, fEMCCut, emcmatchedtracks, nullptr);
    // } else if constexpr (pairtype == PairType::kPCMPHOS) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto phosclusters = std::get<2>(std::tie(args...));
    //   runPairing(collisions, v0photons, phosclusters, v0legs, nullptr, perCollision_pcm, perCollision_phos, fV0PhotonCut, fPHOSCut, nullptr, nullptr);
    // }
    ndf++;
  }
  PROCESS_SWITCH(TaggingPi0, processAnalysis, "process pair analysis", true);

  using FilteredMyCollisionsWithJJMC = soa::Filtered<MyCollisionsWithJJMC>;
  void processAnalysisJJMC(FilteredMyCollisionsWithJJMC const& collisions, Types const&... args)
  {
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == PairType::kPCMDalitzEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runPairing(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut);
    }
    // else if constexpr (pairtype == PairType::kPCMEMC) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto emcclusters = std::get<2>(std::tie(args...));
    //   auto emcmatchedtracks = std::get<3>(std::tie(args...));
    //   runPairing(collisions, v0photons, emcclusters, v0legs, nullptr, perCollision_pcm, perCollision_emc, fV0PhotonCut, fEMCCut, emcmatchedtracks, nullptr);
    // } else if constexpr (pairtype == PairType::kPCMPHOS) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto phosclusters = std::get<2>(std::tie(args...));
    //   runPairing(collisions, v0photons, phosclusters, v0legs, nullptr, perCollision_pcm, perCollision_phos, fV0PhotonCut, fPHOSCut, nullptr, nullptr);
    // }
    ndf++;
  }
  PROCESS_SWITCH(TaggingPi0, processAnalysisJJMC, "process pair analysis", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(TaggingPi0, processDummy, "Dummy function", false);
};
#endif // PWGEM_PHOTONMESON_CORE_TAGGINGPI0_H_
