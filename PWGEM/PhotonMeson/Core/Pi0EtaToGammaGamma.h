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
// This code loops over photons and makes pairs for neutral mesons analyses.
//    Please write to: daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMA_H_
#define PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMA_H_

#include <cstring>
#include <iterator>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <utility>

#include "TString.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PairCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EMTrack.h"
#include "PWGEM/PhotonMeson/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/NMHistograms.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::photonpair;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyPrimaryElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

using MyPrimaryMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMPrimaryMuonsPrefilterBit>;
using MyPrimaryMuon = MyPrimaryMuons::iterator;

using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMCEMEventIds>;
using MyEMCCluster = MyEMCClusters::iterator;

using MyPHOSClusters = soa::Join<aod::PHOSClusters, aod::PHOSEMEventIds>;
using MyPHOSCluster = MyEMCClusters::iterator;

template <PairType pairtype, typename... Types>
struct Pi0EtaToGammaGamma {
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<bool> cfgDoFlow{"cfgDoFlow", false, "flag to analyze vn"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {VARIABLE_WIDTH, -M_PI / 2, -M_PI / 4, 0.0f, +M_PI / 4, +M_PI / 2}, "Mixing bins - event plane angle"};

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
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<bool> cfg_require_v0_on_wwire_ib{"cfg_require_v0_on_wwire_ib", false, "flag to select V0s on W wires ITSib"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.9, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.99, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
    Configurable<bool> cfg_require_v0_with_correct_xz{"cfg_require_v0_with_correct_xz", true, "flag to select V0s with correct xz"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};

    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 10, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 20, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
  } pcmcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.5, "max mass"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.9, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3]"};
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

    // CCDB configuration for PID ML
    Configurable<std::string> BDTLocalPathGamma{"BDTLocalPathGamma", "pid_ml_xgboost.onnx", "Path to the local .onnx file"};

    Configurable<std::string> BDTPathCCDB{"BDTPathCCDB", "Users/d/dsekihat/pwgem/pidml/", "Path on CCDB"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dileptoncuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccut_group";
    Configurable<bool> requireCaloReadout{"requireCaloReadout", true, "Require calorimeters readout when analyzing EMCal/PHOS"};
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

  void init(InitContext&)
  {
    zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
    zvtx_bin_edges.erase(zvtx_bin_edges.begin());

    cent_bin_edges = std::vector<float>(ConfCentBins.value.begin(), ConfCentBins.value.end());
    cent_bin_edges.erase(cent_bin_edges.begin());

    ep_bin_edges = std::vector<float>(ConfEPBins.value.begin(), ConfEPBins.value.end());
    ep_bin_edges.erase(ep_bin_edges.begin());

    emh1 = new o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>(ndepth);
    emh2 = new o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>(ndepth);

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry, cfgDoFlow);
    if constexpr (pairtype == PairType::kPCMDalitzEE) {
      o2::aod::pwgem::photonmeson::utils::nmhistogram::addNMHistograms(&fRegistry, cfgDoFlow, false, "ee#gamma");
    } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
      o2::aod::pwgem::photonmeson::utils::nmhistogram::addNMHistograms(&fRegistry, cfgDoFlow, false, "#mu#mu#gamma");
    } else {
      o2::aod::pwgem::photonmeson::utils::nmhistogram::addNMHistograms(&fRegistry, cfgDoFlow, false, "#gamma#gamma");
    }
    DefineEMEventCut();
    DefinePCMCut();
    DefineDileptonCut();
    DefineEMCCut();
    DefinePHOSCut();

    if constexpr (pairtype == kEMCEMC) {
      fRegistry.addClone("Pair/same/", "Pair/rotation/");
    }
  }

  ~Pi0EtaToGammaGamma()
  {
    delete emh1;
    emh1 = 0x0;
    delete emh2;
    emh2 = 0x0;

    used_photonIds.clear();
    used_photonIds.shrink_to_fit();
    used_dileptonIds.clear();
    used_dileptonIds.shrink_to_fit();
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
    fV0PhotonCut.SetTrackPtRange(pcmcuts.cfg_min_pt_v0 * 0.4, 1e+10f);
    fV0PhotonCut.SetTrackEtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
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

  void DefineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mass, dileptoncuts.cfg_max_mass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.SetPairDCARange(dileptoncuts.cfg_min_pair_dca3d, dileptoncuts.cfg_max_pair_dca3d); // in sigma
    fDileptonCut.ApplyPhiV(dileptoncuts.cfg_apply_phiv);
    fDileptonCut.ApplyPrefilter(dileptoncuts.cfg_apply_pf);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    // for track
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_max_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMeanClusterSizeITSob(0, 16);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaMuRange(dileptoncuts.cfg_min_TPCNsigmaMu, dileptoncuts.cfg_max_TPCNsigmaMu);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTPCNsigmaKaRange(dileptoncuts.cfg_min_TPCNsigmaKa, dileptoncuts.cfg_max_TPCNsigmaKa);
    fDileptonCut.SetTPCNsigmaPrRange(dileptoncuts.cfg_min_TPCNsigmaPr, dileptoncuts.cfg_max_TPCNsigmaPr);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);

    if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      o2::ml::OnnxModel* eid_bdt = new o2::ml::OnnxModel();
      if (dileptoncuts.loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        std::map<std::string, std::string> metadata;
        bool retrieveSuccessGamma = ccdbApi.retrieveBlob(dileptoncuts.BDTPathCCDB.value, ".", metadata, dileptoncuts.timestampCCDB.value, false, dileptoncuts.BDTLocalPathGamma.value);
        if (retrieveSuccessGamma) {
          eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Gamma model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      } else {
        eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
      }

      fDileptonCut.SetPIDModel(eid_bdt);
    } // end of PID ML
  }

  void DefineEMCCut()
  {
    const float a = emccuts.EMC_TM_Eta->at(0);
    const float b = emccuts.EMC_TM_Eta->at(1);
    const float c = emccuts.EMC_TM_Eta->at(2);

    const float d = emccuts.EMC_TM_Phi->at(0);
    const float e = emccuts.EMC_TM_Phi->at(1);
    const float f = emccuts.EMC_TM_Phi->at(2);
    LOGF(info, "EMCal track matching parameters : a = %f, b = %f, c = %f, d = %f, e = %f, f = %f", a, b, c, d, e, f);

    fEMCCut.SetMinE(emccuts.EMC_minE);
    fEMCCut.SetMinNCell(emccuts.EMC_minNCell);
    fEMCCut.SetM02Range(emccuts.EMC_minM02, emccuts.EMC_maxM02);
    fEMCCut.SetTimeRange(emccuts.EMC_minTime, emccuts.EMC_maxTime);

    fEMCCut.SetTrackMatchingEta([&a, &b, &c](float pT) { return a + pow(pT + b, c); });
    fEMCCut.SetTrackMatchingPhi([&d, &e, &f](float pT) { return d + pow(pT + e, f); });

    fEMCCut.SetMinEoverP(emccuts.EMC_Eoverp);
    fEMCCut.SetUseExoticCut(emccuts.EMC_UseExoticCut);
  }

  void DefinePHOSCut()
  {
    fPHOSCut.SetEnergyRange(phoscuts.cfg_min_Ecluster, 1e+10);
  }

  /// \brief Calculate background (using rotation background method only for EMCal!)
  template <typename TPhotons>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, EMCPhotonCut const& cut, aod::SkimEMCMTs const&)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < 3) {
      return;
    }
    const float rotationAngle = M_PI / 2.0; // rotaion angle 90 degree
    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationAngle);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    for (auto& photon : photons_coll) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      if (!cut.template IsSelected<aod::SkimEMCMTs>(photon)) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector photon3(photon.pt(), photon.eta(), photon.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
      ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;

      float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));
      float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));

      if (openingAngle1 > emccuts.minOpenAngle && abs(mother1.Rapidity()) < maxY) {
        fRegistry.fill(HIST("Pair/rotation/hMggPt"), mother1.M(), mother1.Pt());
      }
      if (openingAngle2 > emccuts.minOpenAngle && abs(mother2.Rapidity()) < maxY) {
        fRegistry.fill(HIST("Pair/rotation/hMggPt"), mother2.M(), mother2.Pt());
      }
    }
  }

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;
  Preslice<MyEMCClusters> perCollision_emc = aod::emccluster::emeventId;
  Preslice<MyPHOSClusters> perCollision_phos = aod::phoscluster::emeventId;

  Preslice<MyPrimaryElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Partition<MyPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);
  Partition<MyPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);

  Preslice<MyPrimaryMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Partition<MyPrimaryMuons> muons_pos = o2::aod::emprimarymuon::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaMu) < o2::aod::pidtpc::tpcNSigmaMu&& o2::aod::pidtpc::tpcNSigmaMu < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaMu);
  Partition<MyPrimaryMuons> muons_neg = o2::aod::emprimarymuon::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaMu) < o2::aod::pidtpc::tpcNSigmaMu && o2::aod::pidtpc::tpcNSigmaMu < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaMu);

  o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>* emh1 = nullptr;
  o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>* emh2 = nullptr;
  std::vector<std::pair<int, int>> used_photonIds;              // <ndf, trackId>
  std::vector<std::tuple<int, int, int, int>> used_dileptonIds; // <ndf, trackId>

  template <typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2, typename TTracksMatchedWithEMC, typename TTracksMatchedWithPHOS>
  void runPairing(TCollisions const& collisions,
                  TPhotons1 const& photons1, TPhotons2 const& photons2,
                  TSubInfos1 const& /*subinfos1*/, TSubInfos2 const& /*subinfos2*/,
                  TPreslice1 const& perCollision1, TPreslice2 const& perCollision2,
                  TCut1 const& cut1, TCut2 const& cut2,
                  TTracksMatchedWithEMC const& tracks_emc, TTracksMatchedWithPHOS const& /*tracks_phos*/)
  {
    for (auto& collision : collisions) {
      int ndiphoton = 0;
      if ((pairtype == PairType::kPHOSPHOS || pairtype == PairType::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if ((pairtype == PairType::kEMCEMC || pairtype == PairType::kPCMEMC) && ((!collision.alias_bit(triggerAliases::kTVXinEMC) && emccuts.requireCaloReadout) || collision.ncollsPerBC() != 1)) {
        continue;
      }

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, cfgDoFlow);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, cfgDoFlow);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

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

      float ep2 = collision.ep2ft0c();
      int epbin = lower_bound(ep_bin_edges.begin(), ep_bin_edges.end(), ep2) - ep_bin_edges.begin() - 1;
      if (epbin < 0) {
        epbin = 0;
      } else if (static_cast<int>(ep_bin_edges.size()) - 2 < epbin) {
        epbin = static_cast<int>(ep_bin_edges.size()) - 2;
      }

      std::tuple<int, int, int> key_bin = std::make_tuple(zbin, centbin, epbin);
      std::pair<int, int64_t> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) { // same kinds pairing
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (abs(v12.Rapidity()) > maxY) {
            continue;
          }
          o2::aod::pwgem::photonmeson::utils::nmhistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, v12, cfgDoFlow);

          if constexpr (pairtype == PairType::kEMCEMC) {
            RotationBackground<MyEMCClusters>(v12, v1, v2, photons2_per_collision, g1.globalIndex(), g2.globalIndex(), cut1, tracks_emc);
          }

          std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
          std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, g2.globalIndex());

          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
            emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id1);
          }
          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id2) == used_photonIds.end()) {
            emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g2.globalIndex(), g2.pt(), g2.eta(), g2.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id2);
          }
          ndiphoton++;
        } // end of pairing loop
      } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

        for (auto& g1 : photons1_per_collision) {
          if (!cut1.template IsSelected<TSubInfos1>(g1)) {
            continue;
          }
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

          for (auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {

            if (pos2.trackId() == ele2.trackId()) { // this is protection against pairing identical 2 tracks.
              continue;
            }
            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
              continue;
            }

            if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) {
              if (!cut2.template IsSelectedTrack<true>(pos2, collision) || !cut2.template IsSelectedTrack<true>(ele2, collision)) {
                continue;
              }
            } else { // cut-based
              if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
                continue;
              }
            }

            if (!cut2.IsSelectedPair(pos2, ele2, collision.bz())) {
              continue;
            }

            ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
            ROOT::Math::PtEtaPhiMVector veeg = v_gamma + v_pos + v_ele;
            if (abs(veeg.Rapidity()) > maxY) {
              continue;
            }
            o2::aod::pwgem::photonmeson::utils::nmhistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, veeg, cfgDoFlow);
            float dca_pos_3d = pos2.dca3DinSigma();
            float dca_ele_3d = ele2.dca3DinSigma();
            float dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
            std::tuple<int, int, int, int> tuple_tmp_id2 = std::make_tuple(ndf, collision.globalIndex(), pos2.trackId(), ele2.trackId());
            if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
              emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
              used_photonIds.emplace_back(pair_tmp_id1);
            }
            if (std::find(used_dileptonIds.begin(), used_dileptonIds.end(), tuple_tmp_id2) == used_dileptonIds.end()) {
              emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), -1, v_ee.Pt(), v_ee.Eta(), v_ee.Phi(), v_ee.M(), 0, dca_ee_3d));
              used_dileptonIds.emplace_back(tuple_tmp_id2);
            }
            ndiphoton++;
          } // end of dielectron loop
        }   // end of g1 loop
      } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto muons_pos_per_collision = muons_pos->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
        auto muons_neg_per_collision = muons_neg->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);

        for (auto& g1 : photons1_per_collision) {
          if (!cut1.template IsSelected<TSubInfos1>(g1)) {
            continue;
          }
          auto pos1 = g1.template posTrack_as<TSubInfos1>();
          auto ele1 = g1.template negTrack_as<TSubInfos1>();
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

          for (auto& [muplus, muminus] : combinations(CombinationsFullIndexPolicy(muons_pos_per_collision, muons_neg_per_collision))) {

            if (muplus.trackId() == muminus.trackId()) { // this is protection against pairing identical 2 tracks.
              continue;
            }
            if (pos1.trackId() == muplus.trackId() || ele1.trackId() == muminus.trackId()) {
              continue;
            }

            if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) {
              if (!cut2.template IsSelectedTrack<true>(muplus, collision) || !cut2.template IsSelectedTrack<true>(muminus, collision)) {
                continue;
              }
            } else { // cut-based
              if (!cut2.template IsSelectedTrack<false>(muplus, collision) || !cut2.template IsSelectedTrack<false>(muminus, collision)) {
                continue;
              }
            }

            if (!cut2.IsSelectedPair(muplus, muminus, collision.bz())) {
              continue;
            }

            ROOT::Math::PtEtaPhiMVector v_pos(muplus.pt(), muplus.eta(), muplus.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele(muminus.pt(), muminus.eta(), muminus.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
            ROOT::Math::PtEtaPhiMVector veeg = v_gamma + v_pos + v_ele;
            if (abs(veeg.Rapidity()) > maxY) {
              continue;
            }
            o2::aod::pwgem::photonmeson::utils::nmhistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, veeg, cfgDoFlow);
            float dca_pos_3d = muplus.dca3DinSigma();
            float dca_ele_3d = muminus.dca3DinSigma();
            float dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
            std::tuple<int, int, int, int> tuple_tmp_id2 = std::make_tuple(ndf, collision.globalIndex(), muplus.trackId(), muminus.trackId());
            if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
              emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
              used_photonIds.emplace_back(pair_tmp_id1);
            }
            if (std::find(used_dileptonIds.begin(), used_dileptonIds.end(), tuple_tmp_id2) == used_dileptonIds.end()) {
              emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), -1, v_ee.Pt(), v_ee.Eta(), v_ee.Phi(), v_ee.M(), 0, dca_ee_3d));
              used_dileptonIds.emplace_back(tuple_tmp_id2);
            }
            ndiphoton++;
          } // end of dimuon loop
        }   // end of g1 loop

      } else { // PCM-EMC, PCM-PHOS. Nightmare. don't run these pairs.
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (abs(v12.Rapidity()) > maxY) {
            continue;
          }
          o2::aod::pwgem::photonmeson::utils::nmhistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
          std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
          std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, g2.globalIndex());

          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
            emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id1);
          }
          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id2) == used_photonIds.end()) {
            emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g2.globalIndex(), g2.pt(), g2.eta(), g2.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id2);
          }
          ndiphoton++;
        } // end of pairing loop
      }   // end of pairing in same event

      // event mixing
      if (!cfgDoMix || !(ndiphoton > 0)) {
        continue;
      }

      // make a vector of selected photons in this collision.
      auto selected_photons1_in_this_event = emh1->GetTracksPerCollision(key_df_collision);
      auto selected_photons2_in_this_event = emh2->GetTracksPerCollision(key_df_collision);

      auto collisionIds1_in_mixing_pool = emh1->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIds2_in_mixing_pool = emh2->GetCollisionIdsFromEventPool(key_bin);

      if constexpr (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) { // same kinds pairing
        selected_photons1_in_this_event.insert(selected_photons1_in_this_event.end(), selected_photons2_in_this_event.begin(), selected_photons2_in_this_event.end());

        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);
          auto photons2_from_event_pool = emh2->GetTracksPerCollision(mix_dfId_collisionId);
          photons1_from_event_pool.insert(photons1_from_event_pool.end(), photons2_from_event_pool.begin(), photons2_from_event_pool.end());
          // LOGF(info, "Do event mixing: current event (%d, %d), ngamma = %d | event pool (%d, %d), ngamma = %d", ndf, collision.globalIndex(), selected_photons1_in_this_event.size(), mix_dfId, mix_collisionId, photons1_from_event_pool.size());

          for (auto& g1 : selected_photons1_in_this_event) {
            for (auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              o2::aod::pwgem::photonmeson::utils::nmhistogram::fillPairInfo<1, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
            }
          }
        } // end of loop over mixed event pool

      } else { //[photon1 from event1, photon2 from event2] and [photon1 from event2, photon2 from event1]

        for (auto& mix_dfId_collisionId : collisionIds2_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto photons2_from_event_pool = emh2->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), ngamma = %d | event pool (%d, %d), nll = %d", ndf, collision.globalIndex(), selected_photons1_in_this_event.size(), mix_dfId, mix_collisionId, photons2_from_event_pool.size());

          for (auto& g1 : selected_photons1_in_this_event) {
            for (auto& g2 : photons2_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v2.SetM(g2.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              o2::aod::pwgem::photonmeson::utils::nmhistogram::fillPairInfo<1, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
            }
          }
        } // end of loop over mixed event pool

        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), nll = %d | event pool (%d, %d), ngamma = %d", ndf, collision.globalIndex(), selected_photons2_in_this_event.size(), mix_dfId, mix_collisionId, photons1_from_event_pool.size());

          for (auto& g1 : selected_photons2_in_this_event) {
            for (auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v1.SetM(g1.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              o2::aod::pwgem::photonmeson::utils::nmhistogram::fillPairInfo<1, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
            }
          }
        } // end of loop over mixed event pool
      }

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh2->AddCollisionIdAtLast(key_bin, key_df_collision);
      }

    } // end of collision loop
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  int ndf = 0;
  void processAnalysis(FilteredMyCollisions const& collisions, Types const&... args)
  {
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == PairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      runPairing(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut, nullptr, nullptr);
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runPairing(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut, nullptr, nullptr);
    } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimarymuons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runPairing(collisions, v0photons, emprimarymuons, v0legs, emprimarymuons, perCollision_pcm, perCollision_muon, fV0PhotonCut, fDileptonCut, nullptr, nullptr);
    } else if constexpr (pairtype == PairType::kEMCEMC) {
      auto emcclusters = std::get<0>(std::tie(args...));
      auto emcmatchedtracks = std::get<1>(std::tie(args...));
      runPairing(collisions, emcclusters, emcclusters, nullptr, nullptr, perCollision_emc, perCollision_emc, fEMCCut, fEMCCut, emcmatchedtracks, nullptr);
    } else if constexpr (pairtype == PairType::kPHOSPHOS) {
      auto phosclusters = std::get<0>(std::tie(args...));
      runPairing(collisions, phosclusters, phosclusters, nullptr, nullptr, perCollision_phos, perCollision_phos, fPHOSCut, fPHOSCut, nullptr, nullptr);
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
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processAnalysis, "process pair analysis", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processDummy, "Dummy function", true);
};
#endif // PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMA_H_
