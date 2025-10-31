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
/// \file DiphotonHadronMPC.h
/// \brief This code is to analyze diphoton-hadron correlation. Keep in mind that cumulant method does not require event mixing.
///
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_DIPHOTONHADRONMPC_H_
#define PWGEM_PHOTONMESON_CORE_DIPHOTONHADRONMPC_H_

#include "PWGEM/Dilepton/Core/EMTrackCut.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/NMHistograms.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TString.h"

#include <algorithm>
#include <array>
#include <cstring>
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
using namespace o2::aod::pwgem::photonmeson::photonpair;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithSWT = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMSWTriggerBits>;
using MyCollisionWithSWT = MyCollisionsWithSWT::iterator;

using MyV0Photons = soa::Filtered<soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsKFPrefilterBitDerived>>;
using MyV0Photon = MyV0Photons::iterator;

using MyPrimaryElectrons = soa::Filtered<soa::Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBitDerived>>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

template <PairType pairtype, typename... Types>
struct DiphotonHadronMPC {

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<std::string> cfg_swt_name{"cfg_swt_name", "fHighTrackMult", "desired software trigger name"}; // 1 trigger name per 1 task. fHighTrackMult, fHighFt0Mult

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgOccupancyEstimator{"cfgOccupancyEstimator", 0, "FT0C:0, Track:1"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.8, "maximum rapidity for diphoton"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth_photon{"ndepth_photon", 100, "depth for event mixing between photon-photon"};
  Configurable<int> ndepth_hadron{"ndepth_hadron", 2, "depth for event mixing between hadron-hadron"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 0.1, 1, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  ConfigurableAxis ConfMggBins{"ConfMggBins", {200, 0.0, 0.8}, "mgg bins for output histograms"};
  ConfigurableAxis ConfPtggBins{"ConfPtggBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTgg bins for output histograms"};

  ConfigurableAxis ConfPtHadronBins{"ConfPtHadronBins", {VARIABLE_WIDTH, 0.00, 0.15, 0.2, 0.3, 0.4, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00}, "pT,h bins for output histograms"};
  ConfigurableAxis ConfDEtaBins{"ConfDEtaBins", {120, -6, 6}, "deta bins for output histograms"};
  Configurable<int> cfgNbinsDPhi{"cfgNbinsDPhi", 36, "nbins in dphi for output histograms"};
  // Configurable<int> cfgNbinsCosNDPhi{"cfgNbinsCosNDPhi", 100, "nbins in cos(n(dphi)) for output histograms"};
  // Configurable<int> cfgNmod{"cfgNmod", 2, "n-th harmonics"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    // for RCT
    Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
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
    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to V0"};

    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};
  } pcmcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.04, "max mass"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 2.0, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.1, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.1, "max dca Z for single track in cm"};
    Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1e+10, "max DCA 3D in sigma"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to electron"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTOFif), "pid scheme [kTOFif : 0, kTPConly : 1]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -0.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +0.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
  } dileptoncuts;

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
    Configurable<uint16_t> cfg_track_bits{"cfg_track_bits", 645, "required track bits"}; // default:645, loose:0, tight:778
  } trackcuts;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> zvtx_bin_edges;
  std::vector<float> cent_bin_edges;
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

    LOGF(info, "cfgOccupancyEstimator = %d", cfgOccupancyEstimator.value);
    occ_bin_edges = std::vector<float>(ConfOccupancyBins.value.begin(), ConfOccupancyBins.value.end());
    occ_bin_edges.erase(occ_bin_edges.begin());

    emh1 = new MyEMH(ndepth_photon);
    emh2 = new MyEMH(ndepth_photon);
    emh_diphoton = new MyEMH_track(ndepth_photon);
    emh_ref = new MyEMH_track(ndepth_hadron);

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    addHistograms();

    DefineEMEventCut();
    DefineEMTrackCut();
    DefinePCMCut();
    DefineDileptonCut();

    fRegistry.add("Diphoton/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", kTH1D, {{10001, -0.5, 10000.5}}, true);
    if (doprocessTriggerAnalysis) {
      fRegistry.add("Event/hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
    }

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    rctChecker.init(eventcuts.cfgRCTLabel.value, eventcuts.cfgCheckZDC.value, eventcuts.cfgTreatLimitedAcceptanceAsBad.value);
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
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    }
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
      // LOGF(info, "total inspected TVX events = %d in run number %d", collision.nInspectedTVX(), collision.runNumber());
      // fRegistry.fill(HIST("Event/hNInspectedTVX"), collision.runNumber(), collision.nInspectedTVX());
    }
  }

  ~DiphotonHadronMPC()
  {
    delete emh1;
    emh1 = 0x0;
    delete emh2;
    emh2 = 0x0;
    delete emh_diphoton;
    emh_diphoton = 0x0;
    delete emh_ref;
    emh_ref = 0x0;

    used_photonIds_per_col.clear();
    used_photonIds_per_col.shrink_to_fit();
    used_dileptonIds_per_col.clear();
    used_dileptonIds_per_col.shrink_to_fit();
    used_diphotonIds_per_col.clear();
    used_diphotonIds_per_col.shrink_to_fit();

    map_mixed_eventId_to_globalBC.clear();
  }

  void addHistograms()
  {
    std::string mass_axis_title = "m_{#gamma#gamma} (GeV/c^{2})";
    std::string pair_pt_axis_title = "p_{T,#gamma#gamma} (GeV/c)";
    std::string deta_axis_title = "#Delta#eta = #eta_{#gamma#gamma} - #eta_{h}";
    std::string dphi_axis_title = "#Delta#varphi = #varphi_{#gamma#gamma} - #varphi_{h} (rad.)";
    // std::string cosndphi_axis_title = std::format("cos({0:d}(#varphi_{{#gamma#gamma}} - #varphi_{{h}}))", cfgNmod.value);

    if constexpr (pairtype == PairType::kPCMDalitzEE) {
      mass_axis_title = "m_{ee#gamma} (GeV/c^{2})";
      pair_pt_axis_title = "p_{T,ee#gamma} (GeV/c)";
      deta_axis_title = "#Delta#eta = #eta_{ee#gamma} - #eta_{h}";
      dphi_axis_title = "#Delta#varphi = #varphi_{ee#gamma} - #varphi_{h} (rad.)";
      // cosndphi_axis_title = std::format("cos({0:d}(#varphi_{{ee#gamma}} - #varphi_{{h}}))", cfgNmod.value);
    }

    // diphoton info
    const AxisSpec axis_mass{ConfMggBins, mass_axis_title};
    const AxisSpec axis_pt{ConfPtggBins, pair_pt_axis_title};

    // diphoton-hadron info
    const AxisSpec axis_deta{ConfDEtaBins, deta_axis_title};
    const AxisSpec axis_dphi{cfgNbinsDPhi, -M_PI / 2, +3 * M_PI / 2, dphi_axis_title};

    const AxisSpec axis_pt_hadron{ConfPtHadronBins, "p_{T,h} (GeV/c)"};
    const AxisSpec axis_eta_hadron{40, -2, +2, "#eta_{h}"};
    const AxisSpec axis_phi_hadron{36, 0, 2 * M_PI, "#varphi_{h} (rad.)"};

    fRegistry.add("Hadron/hs", "hadron", kTHnSparseD, {axis_pt_hadron, axis_eta_hadron, axis_phi_hadron}, false);
    fRegistry.add("Hadron/hTrackBit", "track bit", kTH1D, {{65536, -0.5, 65535.5}}, false);

    fRegistry.add("Diphoton/same/hs", "diphoton", kTHnSparseD, {axis_mass, axis_pt}, true);
    fRegistry.addClone("Diphoton/same/", "Diphoton/mix/");

    fRegistry.add("DiphotonHadron/same/hs", "diphoton-hadron 2PC", kTHnSparseD, {axis_mass, axis_pt, axis_deta, axis_dphi}, true);
    fRegistry.addClone("DiphotonHadron/same/", "DiphotonHadron/mix/");

    // hadron-hadron
    const AxisSpec axis_deta_hh{60, -3, +3, "#Delta#eta = #eta_{h}^{ref1} - #eta_{h}^{ref2}"};
    const AxisSpec axis_dphi_hh{90, -M_PI / 2, +3 * M_PI / 2, "#Delta#varphi = #varphi_{h}^{ref1} - #varphi_{h}^{ref2} (rad.)"};
    // const AxisSpec axis_cosndphi_hh{cfgNbinsCosNDPhi, -1, +1, std::format("cos({0:d}(#varphi_{{h}}^{{ref1}} - #varphi_{{h}}^{{ref2}}))", cfgNmod.value)};
    fRegistry.add("HadronHadron/same/hDEtaDPhi", "hadron-hadron 2PC", kTH2D, {axis_dphi_hh, axis_deta_hh}, true);
    fRegistry.addClone("HadronHadron/same/", "HadronHadron/mix/");
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(eventcuts.cfgZvtxMin, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
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
    fDileptonCut.IncludeITSsa(false, 0.15);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);
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

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;

  Preslice<MyPrimaryElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Partition<MyPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);
  Partition<MyPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);

  using RefTracks = soa::Join<aod::EMPrimaryTracks, aod::EMPrimaryTrackEMEventIds>;
  using RefTrack = RefTracks::iterator;
  Preslice<RefTracks> perCollision_track = aod::emprimarytrack::emeventId;
  Filter refTrackFilter = trackcuts.cfg_min_pt_track < 1 / nabs(o2::aod::emprimarytrack::signed1Pt) && 1 / nabs(o2::aod::emprimarytrack::signed1Pt) < trackcuts.cfg_max_pt_track && trackcuts.cfg_min_eta_track < o2::aod::emprimarytrack::eta && o2::aod::emprimarytrack::eta < trackcuts.cfg_max_eta_track;
  using FilteredRefTracks = soa::Filtered<RefTracks>;
  using FilteredRefTrack = FilteredRefTracks::iterator;

  using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrack>;
  MyEMH* emh1 = nullptr;
  MyEMH* emh2 = nullptr;
  using MyEMH_track = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrack>; // for charged track
  MyEMH_track* emh_diphoton = nullptr;
  MyEMH_track* emh_ref = nullptr;

  std::vector<int> used_photonIds_per_col;                         // <ndf, trackId>
  std::vector<std::pair<int, int>> used_dileptonIds_per_col;       // <ndf, trackId>
  std::vector<std::tuple<int, int, int>> used_diphotonIds_per_col; // <ndf, trackId>
  std::map<std::pair<int, int>, uint64_t> map_mixed_eventId_to_globalBC;

  template <bool isTriggerAnalysis, typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2, typename TRefTracks>
  void run2PC(TCollisions const& collisions,
              TPhotons1 const& photons1, TPhotons2 const& photons2, TSubInfos1 const&, TSubInfos2 const&, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCut1 const& cut1, TCut2 const& cut2,
              TRefTracks const& refTracks)
  {
    for (const auto& collision : collisions) {
      initCCDB<isTriggerAnalysis>(collision);
      int ndiphoton = 0;

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (!collision.swtalias_bit(o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value))) {
          continue;
        }
      }

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, 1.f);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (eventcuts.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, 1.f);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0, 1.f); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0, 1.f);  // accepted

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

      // LOGF(info, "collision.globalIndex() = %d, collision.posZ() = %f, centrality = %f, ep2 = %f, collision.trackOccupancyInTimeRange() = %d, zbin = %d, centbin = %d, epbin = %d, occbin = %d", collision.globalIndex(), collision.posZ(), centrality, ep2, collision.trackOccupancyInTimeRange(), zbin, centbin, epbin, occbin);

      auto refTracks_per_collision = refTracks.sliceBy(perCollision_track, collision.globalIndex());

      std::tuple<int, int, int, int> key_bin = std::make_tuple(zbin, centbin, epbin, occbin);
      std::pair<int, int64_t> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == PairType::kPCMPCM) { // same kinds pairing
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (std::fabs(v12.Rapidity()) > maxY) {
            continue;
          }
          fRegistry.fill(HIST("Diphoton/same/hs"), v12.M(), v12.Pt());
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          auto pos2 = g2.template posTrack_as<TSubInfos2>();
          auto ele2 = g2.template negTrack_as<TSubInfos2>();

          int npair = 0; // the number of diphoton-h pairs
          for (const auto& track : refTracks_per_collision) {
            if (pos1.trackId() == track.trackId() || ele1.trackId() == track.trackId()) {
              continue;
            }
            if (pos2.trackId() == track.trackId() || ele2.trackId() == track.trackId()) {
              continue;
            }

            if (fEMTrackCut.IsSelected(track)) {
              ROOT::Math::PtEtaPhiMVector v3(track.pt(), track.eta(), track.phi(), 0.139);
              float deta = v12.Eta() - v3.Eta();
              float dphi = v12.Phi() - v3.Phi();
              // o2::math_utils::bringTo02Pi(dphi);
              dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
              fRegistry.fill(HIST("DiphotonHadron/same/hs"), v12.M(), v12.Pt(), deta, dphi);
              npair++;
            }
          } // end of ref track loop

          if (npair > 0) {
            std::tuple<int, int, int> tuple_tmp_diphoton = std::make_tuple(g1.globalIndex(), g2.globalIndex(), -1);
            if (std::find(used_diphotonIds_per_col.begin(), used_diphotonIds_per_col.end(), tuple_tmp_diphoton) == used_diphotonIds_per_col.end()) {
              emh_diphoton->AddTrackToEventPool(key_df_collision, EMTrack(v12.Pt(), v12.Eta(), v12.Phi(), v12.M()));
              used_diphotonIds_per_col.emplace_back(tuple_tmp_diphoton);
            }
          }

          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g1.globalIndex()) == used_photonIds_per_col.end()) {
            emh1->AddTrackToEventPool(key_df_collision, EMTrack(g1.pt(), g1.eta(), g1.phi(), 0));
            used_photonIds_per_col.emplace_back(g1.globalIndex());
          }
          if (std::find(used_photonIds_per_col.begin(), used_photonIds_per_col.end(), g2.globalIndex()) == used_photonIds_per_col.end()) {
            emh1->AddTrackToEventPool(key_df_collision, EMTrack(g2.pt(), g2.eta(), g2.phi(), 0));
            used_photonIds_per_col.emplace_back(g2.globalIndex());
          }
          ndiphoton++;
        } // end of pairing loop
      } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

        for (const auto& g1 : photons1_per_collision) {
          if (!cut1.template IsSelected<TSubInfos1>(g1)) {
            continue;
          }
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

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
            if (std::fabs(veeg.Rapidity()) > maxY) {
              continue;
            }
            fRegistry.fill(HIST("Diphoton/same/hs"), veeg.M(), veeg.Pt());

            int npair = 0; // the number of diphoton-h pairs
            for (const auto& track : refTracks_per_collision) {
              if (pos1.trackId() == track.trackId() || ele1.trackId() == track.trackId()) {
                continue;
              }
              if (pos2.trackId() == track.trackId() || ele2.trackId() == track.trackId()) {
                continue;
              }

              ROOT::Math::PtEtaPhiMVector v3(track.pt(), track.eta(), track.phi(), 0.139);
              float deta = veeg.Eta() - v3.Eta();
              float dphi = veeg.Phi() - v3.Phi();
              // o2::math_utils::bringTo02Pi(dphi);
              dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
              fRegistry.fill(HIST("DiphotonHadron/same/hs"), veeg.M(), veeg.Pt(), deta, dphi);
              npair++;

            } // end of ref track loop

            if (npair > 0) {
              std::tuple<int, int, int> tuple_tmp_diphoton = std::make_tuple(g1.globalIndex(), pos2.trackId(), ele2.trackId());
              if (std::find(used_diphotonIds_per_col.begin(), used_diphotonIds_per_col.end(), tuple_tmp_diphoton) == used_diphotonIds_per_col.end()) {
                emh_diphoton->AddTrackToEventPool(key_df_collision, EMTrack(veeg.Pt(), veeg.Eta(), veeg.Phi(), veeg.M()));
                used_diphotonIds_per_col.emplace_back(tuple_tmp_diphoton);
              }
            }

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
      } // end of pairing in same event

      used_photonIds_per_col.clear();
      used_photonIds_per_col.shrink_to_fit();
      used_dileptonIds_per_col.clear();
      used_dileptonIds_per_col.shrink_to_fit();
      used_diphotonIds_per_col.clear();
      used_diphotonIds_per_col.shrink_to_fit();

      if (ndiphoton > 0) {
        emh_ref->ReserveNTracksPerCollision(key_df_collision, refTracks_per_collision.size());
        for (const auto& track : refTracks_per_collision) {
          if (fEMTrackCut.IsSelected(track)) {
            fRegistry.fill(HIST("Hadron/hs"), track.pt(), track.eta(), track.phi());
            fRegistry.fill(HIST("Hadron/hTrackBit"), track.trackBit());
            emh_ref->AddTrackToEventPool(key_df_collision, EMTrack(track.pt(), track.eta(), track.phi(), 0.139));
          }
        }

        for (const auto& [ref1, ref2] : combinations(CombinationsStrictlyUpperIndexPolicy(refTracks_per_collision, refTracks_per_collision))) {
          if (fEMTrackCut.IsSelected(ref1) && fEMTrackCut.IsSelected(ref2)) {
            float deta = ref1.eta() - ref2.eta();
            float dphi = ref1.phi() - ref2.phi();
            // o2::math_utils::bringTo02Pi(dphi);
            dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
            fRegistry.fill(HIST("HadronHadron/same/hDEtaDPhi"), dphi, deta);
          }
        }
      }

      // event mixing
      if (!cfgDoMix || !(ndiphoton > 0)) {
        continue;
      }

      // make a vector of selected photons in this collision.
      auto selected_photons1_in_this_event = emh1->GetTracksPerCollision(key_df_collision);
      auto selected_photons2_in_this_event = emh2->GetTracksPerCollision(key_df_collision);
      auto selected_refTracks_in_this_event = emh_ref->GetTracksPerCollision(key_df_collision);
      auto selected_diphotons_in_this_event = emh_diphoton->GetTracksPerCollision(key_df_collision);

      auto collisionIds1_in_mixing_pool = emh1->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIds2_in_mixing_pool = emh2->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIdsRef_in_mixing_pool = emh_ref->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIdsDiphoton_in_mixing_pool = emh_diphoton->GetCollisionIdsFromEventPool(key_bin);

      if constexpr (pairtype == PairType::kPCMPCM) { // same kinds pairing
        for (const auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("Diphoton/mix/hDiffBC"), diffBC);
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
              fRegistry.fill(HIST("Diphoton/mix/hs"), v12.M(), v12.Pt(), 1.f);
            }
          }
        } // end of loop over mixed event pool between photon-photon

        for (const auto& mix_dfId_collisionId : collisionIdsRef_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto refTracks_from_event_pool = emh_ref->GetTracksPerCollision(mix_dfId_collisionId);
          for (const auto& trg : selected_diphotons_in_this_event) {
            for (const auto& ref : refTracks_from_event_pool) {
              float deta = trg.eta() - ref.eta();
              float dphi = trg.phi() - ref.phi();
              // o2::math_utils::bringTo02Pi(dphi);
              dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
              fRegistry.fill(HIST("DiphotonHadron/mix/hs"), trg.mass(), trg.pt(), deta, dphi);
            }
          }
        } // end of loop over mixed event pool between diphoton-hadron

      } else { // [photon1 from event1, photon2 from event2] and [photon1 from event2, photon2 from event1]
        for (const auto& mix_dfId_collisionId : collisionIds2_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("Diphoton/mix/hDiffBC"), diffBC);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto photons2_from_event_pool = emh2->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), ngamma = %d | event pool (%d, %d), nll = %d", ndf, collision.globalIndex(), selected_photons1_in_this_event.size(), mix_dfId, mix_collisionId, photons2_from_event_pool.size());

          for (const auto& g1 : selected_photons1_in_this_event) {
            for (const auto& g2 : photons2_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == PairType::kPCMDalitzEE) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v2.SetM(g2.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (std::fabs(v12.Rapidity()) > maxY) {
                continue;
              }
              fRegistry.fill(HIST("Diphoton/mix/hs"), v12.M(), v12.Pt(), 1.f);
            }
          }
        } // end of loop over mixed event pool between photon-photon

        for (const auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          fRegistry.fill(HIST("Diphoton/mix/hDiffBC"), diffBC);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), nll = %d | event pool (%d, %d), ngamma = %d", ndf, collision.globalIndex(), selected_photons2_in_this_event.size(), mix_dfId, mix_collisionId, photons1_from_event_pool.size());

          for (const auto& g1 : selected_photons2_in_this_event) {
            for (const auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == PairType::kPCMDalitzEE) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v1.SetM(g1.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (std::fabs(v12.Rapidity()) > maxY) {
                continue;
              }
              fRegistry.fill(HIST("Diphoton/mix/hs"), v12.M(), v12.Pt(), 1.f);
            }
          }
        } // end of loop over mixed event pool between photon-photon

        for (const auto& mix_dfId_collisionId : collisionIdsRef_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
          uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
          if (diffBC < ndiff_bc_mix) {
            continue;
          }

          auto refTracks_from_event_pool = emh_ref->GetTracksPerCollision(mix_dfId_collisionId);
          for (const auto& trg : selected_diphotons_in_this_event) {
            for (const auto& ref : refTracks_from_event_pool) {
              float deta = trg.eta() - ref.eta();
              float dphi = trg.phi() - ref.phi();
              // o2::math_utils::bringTo02Pi(dphi);
              dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
              fRegistry.fill(HIST("DiphotonHadron/mix/hs"), trg.mass(), trg.pt(), deta, dphi);
            }
          }
        } // end of loop over mixed event pool between diphoton-hadron
      }

      // hadron-hadron mixed event
      for (const auto& mix_dfId_collisionId : collisionIdsRef_in_mixing_pool) {
        int mix_dfId = mix_dfId_collisionId.first;
        int64_t mix_collisionId = mix_dfId_collisionId.second;

        if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
          continue;
        }

        auto globalBC_mix = map_mixed_eventId_to_globalBC[mix_dfId_collisionId];
        uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
        if (diffBC < ndiff_bc_mix) {
          continue;
        }

        auto refTracks_from_event_pool = emh_ref->GetTracksPerCollision(mix_dfId_collisionId);
        for (const auto& ref1 : selected_refTracks_in_this_event) {
          for (const auto& ref2 : refTracks_from_event_pool) {
            float deta = ref1.eta() - ref2.eta();
            float dphi = ref1.phi() - ref2.phi();
            // o2::math_utils::bringTo02Pi(dphi);
            dphi = RecoDecay::constrainAngle(dphi, -M_PI / 2, 1U);
            fRegistry.fill(HIST("HadronHadron/mix/hDEtaDPhi"), dphi, deta);
          }
        }
      } // end of loop over mixed event pool between hadron-hadron

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh2->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_diphoton->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_ref->AddCollisionIdAtLast(key_bin, key_df_collision);
        map_mixed_eventId_to_globalBC[key_df_collision] = collision.globalBC();
      }

    } // end of collision loop
  }

  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Filter prefilter_pcm = ifnode(pcmcuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::v0photonkf::pfbderived == static_cast<uint16_t>(0), true);
  Filter prefilter_primaryelectron = ifnode(dileptoncuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::emprimaryelectron::pfbderived == static_cast<uint16_t>(0), true);

  int ndf = 0;
  void processAnalysis(FilteredMyCollisions const& collisions, FilteredRefTracks const& refTracks, Types const&... args)
  {
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == PairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      run2PC<false>(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut, refTracks);
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      run2PC<false>(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut, refTracks);
    }
    ndf++;
  }
  PROCESS_SWITCH(DiphotonHadronMPC, processAnalysis, "process pair analysis", true);

  // using FilteredMyCollisionsWithSWT = soa::Filtered<MyCollisionsWithSWT>;
  void processTriggerAnalysis(MyCollisionsWithSWT const& collisions, FilteredRefTracks const& refTracks, Types const&... args)
  {
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == PairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      run2PC<true>(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut, refTracks);
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      run2PC<true>(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut, refTracks);
    }
    ndf++;
  }
  PROCESS_SWITCH(DiphotonHadronMPC, processTriggerAnalysis, "process pair analysis with software trigger", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(DiphotonHadronMPC, processDummy, "Dummy function", false);
};
#endif // PWGEM_PHOTONMESON_CORE_DIPHOTONHADRONMPC_H_
