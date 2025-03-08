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

#ifndef PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMAMC_H_
#define PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMAMC_H_

#include <string>
#include <map>
#include <vector>

#include "TF1.h"
#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/NMHistograms.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::photonpair;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMMCEventLabels, aod::EMEventsWeight>;
using MyCollision = MyCollisions::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::BinnedGenPts>;
using MyMCCollision = MyMCCollisions::iterator;

using MyV0Photons = soa::Filtered<soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsKFPrefilterBitDerived>>;
using MyV0Photon = MyV0Photons::iterator;

using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMEMCClusterMCLabels, aod::EMCEMEventIds>;
using MyEMCCluster = MyEMCClusters::iterator;

using MyPHOSClusters = soa::Join<aod::PHOSClusters, aod::PHOSEMEventIds>;
using MyPHOSCluster = MyEMCClusters::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

using MyMCElectrons = soa::Filtered<soa::Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBitDerived, aod::EMPrimaryElectronMCLabels>>;
using MyMCElectron = MyMCElectrons::iterator;

template <PairType pairtype, typename... Types>
struct Pi0EtaToGammaGammaMC {
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<float> maxY_rec{"maxY_rec", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<std::string> fd_k0s_to_pi0{"fd_k0s_pi0", "1.0", "feed down correction to pi0"};
  Configurable<bool> cfgRequireTrueAssociation{"cfgRequireTrueAssociation", false, "flag to require true mc collision association"};

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
    Configurable<bool> cfg_require_v0_with_correct_xz{"cfg_require_v0_with_correct_xz", true, "flag to select V0s with correct xz"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to V0"};

    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 10, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
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
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.05, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.05, "max dca Z for single track in cm"};
    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to electron"};

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
  } emccuts;

  Configurable<float> maxY_gen{"maxY_gen", 0.9, "maximum rapidity for generated particles"}; // for PCM and dielectron
  Configurable<float> maxRgen{"maxRgen", 90.f, "maximum radius for generated particles"};
  Configurable<float> margin_z_mc{"margin_z_mc", 7.0, "margin for z cut in cm for MC"};

  PHOSPhotonCut fPHOSCut;
  struct : ConfigurableGroup {
    std::string prefix = "phoscut_group";
    Configurable<float> cfg_min_Ecluster{"cfg_min_Ecluster", 0.3, "Minimum cluster energy for PHOS in GeV"};
  } phoscuts;

  TF1* f1fd_k0s_to_pi0;
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};
  static constexpr std::string_view parnames[2] = {"Pi0/", "Eta/"};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  void init(InitContext&)
  {
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    if constexpr (pairtype == PairType::kPCMDalitzEE) {
      o2::aod::pwgem::photonmeson::utils::nmhistogram::addNMHistograms(&fRegistry, true, "ee#gamma");
    } else {
      o2::aod::pwgem::photonmeson::utils::nmhistogram::addNMHistograms(&fRegistry, true, "#gamma#gamma");
    }
    DefineEMEventCut();
    DefinePCMCut();
    DefineDileptonCut();
    DefineEMCCut();
    DefinePHOSCut();

    f1fd_k0s_to_pi0 = new TF1("f1fd_k0s_to_pi0", TString(fd_k0s_to_pi0), 0.f, 100.f);

    fRegistry.add("Event/hNrecPerMCCollision", "Nrec per mc collision;N_{rec} collisions per MC collision", kTH1F, {{21, -0.5f, 20.5f}}, false);

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

  ~Pi0EtaToGammaGammaMC()
  {
    delete f1fd_k0s_to_pi0;
    f1fd_k0s_to_pi0 = 0x0;
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
    fV0PhotonCut.SetTrackPtRange(pcmcuts.cfg_min_pt_v0 * 0.4, 1e+10f);
    fV0PhotonCut.SetTrackEtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfg_max_frac_shared_clusters_tpc);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetIsWithinBeamPipe(pcmcuts.cfg_require_v0_with_correct_xz);
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

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);
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

    fEMCCut = EMCPhotonCut("fEMCCut", "fEMCCut");

    fEMCCut.SetClusterizer(emccuts.clusterDefinition);
    fEMCCut.SetMinE(emccuts.EMC_minE);
    fEMCCut.SetMinNCell(emccuts.EMC_minNCell);
    fEMCCut.SetM02Range(emccuts.EMC_minM02, emccuts.EMC_maxM02);
    fEMCCut.SetTimeRange(emccuts.EMC_minTime, emccuts.EMC_maxTime);

    fEMCCut.SetTrackMatchingEta([a, b, c](float pT) { return a + std::pow(pT + b, c); });
    fEMCCut.SetTrackMatchingPhi([d, e, f](float pT) { return d + std::pow(pT + e, f); });

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

  Preslice<MyMCElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Partition<MyMCElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);
  Partition<MyMCElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);

  // Preslice<MyMCMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  // Partition<MyMCMuons> muons_pos = o2::aod::emprimarymuon::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaMu) < o2::aod::pidtpc::tpcNSigmaMu&& o2::aod::pidtpc::tpcNSigmaMu < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaMu);
  // Partition<MyMCMuons> muons_neg = o2::aod::emprimarymuon::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaMu) < o2::aod::pidtpc::tpcNSigmaMu && o2::aod::pidtpc::tpcNSigmaMu < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaMu);

  template <typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2, typename TMCCollisions, typename TMCParticles>
  void runTruePairing(TCollisions const& collisions,
                      TPhotons1 const& photons1, TPhotons2 const& photons2,
                      TSubInfos1 const& /*subinfos1*/, TSubInfos2 const& /*subinfos2*/,
                      TPreslice1 const& perCollision1, TPreslice2 const& perCollision2,
                      TCut1 const& cut1, TCut2 const& cut2,
                      TMCCollisions const& mccollisions, TMCParticles const& mcparticles)
  {
    for (auto& collision : collisions) {
      initCCDB(collision);

      if ((pairtype == PairType::kPHOSPHOS || pairtype == PairType::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }

      if (eventcuts.onlyKeepWeightedEvents && std::fabs(collision.weight() - 1.) < 1E-10) {
        continue;
      }

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, collision.weight());
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, collision.weight());
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0, collision.weight()); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0, collision.weight());  // accepted

      int photonid1 = -1, photonid2 = -1, pi0id = -1, etaid = -1;
      if constexpr (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) { // same kinds pairing
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }

          if constexpr (pairtype == PairType::kPCMPCM) { // check 2 legs
            auto pos1 = g1.template posTrack_as<TSubInfos1>();
            auto ele1 = g1.template negTrack_as<TSubInfos1>();
            auto pos2 = g2.template posTrack_as<TSubInfos2>();
            auto ele2 = g2.template negTrack_as<TSubInfos2>();

            auto pos1mc = pos1.template emmcparticle_as<TMCParticles>();
            auto ele1mc = ele1.template emmcparticle_as<TMCParticles>();
            auto pos2mc = pos2.template emmcparticle_as<TMCParticles>();
            auto ele2mc = ele2.template emmcparticle_as<TMCParticles>();

            photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
            photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles);
          } else if constexpr (pairtype == PairType::kEMCEMC) {
            auto cluster1mcparticle = mcparticles.iteratorAt(g1.emmcparticleId());
            auto cluster2mcparticle = mcparticles.iteratorAt(g2.emmcparticleId());

            photonid1 = FindMotherInChain(cluster1mcparticle, mcparticles, std::vector<int>{111, 221});
            photonid2 = FindMotherInChain(cluster2mcparticle, mcparticles, std::vector<int>{111, 221});
          } else {
            photonid1 = -1;
            photonid2 = -1;
          }

          if (photonid1 < 0 || photonid2 < 0) {
            continue;
          }
          auto g1mc = mcparticles.iteratorAt(photonid1);
          auto g2mc = mcparticles.iteratorAt(photonid2);

          if constexpr (pairtype == PairType::kPCMPCM) {
            if (!IsConversionPointInAcceptance(g1mc, maxRgen, maxY_gen, margin_z_mc, mcparticles) || !IsConversionPointInAcceptance(g2mc, maxRgen, maxY_gen, margin_z_mc, mcparticles)) {
              continue;
            }
          }

          pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
          etaid = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 221, mcparticles);

          if (g1mc.globalIndex() != g2mc.globalIndex() && pi0id < 0 && etaid < 0) { // for same gamma no pi0/eta will be found, but we still want to fill the FromSameGamma hist
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (std::fabs(v12.Rapidity()) > maxY_rec) {
            continue;
          }

          if (pairtype == PairType::kEMCEMC) {
            float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
            if (openingAngle < emccuts.minOpenAngle) {
              continue;
            }
          }

          if (g1mc.globalIndex() == g2mc.globalIndex()) {
            if (getMotherPDGCode(g1mc, mcparticles) == 111)
              fRegistry.fill(HIST("Pair/Pi0/hs_FromSameGamma"), v12.M(), v12.Pt(), collision.weight());
            else if (getMotherPDGCode(g1mc, mcparticles) == 221)
              fRegistry.fill(HIST("Pair/Eta/hs_FromSameGamma"), v12.M(), v12.Pt(), collision.weight());
            continue;
          }

          if (pi0id > 0) {
            auto pi0mc = mcparticles.iteratorAt(pi0id);
            if (cfgRequireTrueAssociation && (pi0mc.emmceventId() != collision.emmceventId())) {
              continue;
            }
            o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, pi0mc, mcparticles, mccollisions, f1fd_k0s_to_pi0, collision.weight());
          } else if (etaid > 0) {
            auto etamc = mcparticles.iteratorAt(etaid);
            if (cfgRequireTrueAssociation && (etamc.emmceventId() != collision.emmceventId())) {
              continue;
            }
            o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, etamc, mcparticles, mccollisions, f1fd_k0s_to_pi0, collision.weight());
          }
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
          auto pos1mc = pos1.template emmcparticle_as<TMCParticles>();
          auto ele1mc = ele1.template emmcparticle_as<TMCParticles>();
          photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
          if (photonid1 < 0) {
            continue;
          }
          auto g1mc = mcparticles.iteratorAt(photonid1);
          if (!IsConversionPointInAcceptance(g1mc, maxRgen, maxY_gen, margin_z_mc, mcparticles)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.f);

          for (auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) { // ULS
            if (pos2.trackId() == ele2.trackId()) {                                                                                // this is protection against pairing identical 2 tracks.
              continue;
            }
            if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
              continue;
            }

            if (!cut2.template IsSelectedTrack<false>(pos2, collision) || !cut2.template IsSelectedTrack<false>(ele2, collision)) {
              continue;
            }

            if (!cut2.template IsSelectedPair(pos2, ele2, d_bz)) {
              continue;
            }

            auto pos2mc = mcparticles.iteratorAt(pos2.emmcparticleId());
            auto ele2mc = mcparticles.iteratorAt(ele2.emmcparticleId());
            pi0id = FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 111, mcparticles);
            etaid = FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 221, mcparticles);
            if (pi0id < 0 && etaid < 0) {
              continue;
            }
            ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector veeg = v_gamma + v_pos + v_ele;
            if (std::fabs(veeg.Rapidity()) > maxY_rec) {
              continue;
            }
            if (pi0id > 0) {
              auto pi0mc = mcparticles.iteratorAt(pi0id);
              if (cfgRequireTrueAssociation && (pi0mc.emmceventId() != collision.emmceventId())) {
                continue;
              }
              o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, veeg, pi0mc, mcparticles, mccollisions, f1fd_k0s_to_pi0, collision.weight());
            } else if (etaid > 0) {
              auto etamc = mcparticles.iteratorAt(etaid);
              if (cfgRequireTrueAssociation && (etamc.emmceventId() != collision.emmceventId())) {
                continue;
              }
              o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, veeg, etamc, mcparticles, mccollisions, f1fd_k0s_to_pi0, collision.weight());
            }
          } // end of dielectron loop
        } // end of pcm loop
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
          if (std::fabs(v12.Rapidity()) > maxY_rec) {
            continue;
          }
          // if (pi0id > 0) {
          //   auto pi0mc = mcparticles.iteratorAt(pi0id);
          //   o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, pi0mc, mcparticles, mccollisions, f1fd_k0s_to_pi0, collision.weight());
          // } else if (etaid > 0) {
          //   auto etamc = mcparticles.iteratorAt(etaid);
          //   o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, etamc, mcparticles, mccollisions, f1fd_k0s_to_pi0, collision.weight());
          // }
        } // end of pairing loop
      } // end of pairing in same event
    } // end of collision loop
  }

  template <int par_id, typename TBinnedData>
  void fillBinnedData(TBinnedData const& binned_data, const float weight = 1.f)
  {
    int xbin = 0, ybin = 0, zbin = 0;
    auto hPtY = fRegistry.get<TH2>(HIST("Generated/") + HIST(parnames[par_id]) + HIST("hPtY")); // 2D
    auto hPt = fRegistry.get<TH1>(HIST("Generated/") + HIST(parnames[par_id]) + HIST("hPt"));   // 1D

    for (int ibin = 0; ibin < hPtY->GetNcells(); ibin++) {
      int nentry = binned_data[ibin];
      hPtY->GetBinXYZ(ibin, xbin, ybin, zbin);
      float pt = hPtY->GetXaxis()->GetBinCenter(xbin);
      float y = hPtY->GetYaxis()->GetBinCenter(ybin);
      if (y > maxY_gen) {
        continue;
      }

      for (int j = 0; j < nentry; j++) {
        hPtY->Fill(pt, y, weight);
        hPt->Fill(pt, weight);
      }
    }
  }

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  PresliceUnsorted<MyCollisions> rec_perMcCollision = aod::emmceventlabel::emmceventId;

  template <typename TCollisions, typename TMCCollisions, typename TMCParticles>
  void runGenInfo(TCollisions const& collisions, TMCCollisions const& mccollisions, TMCParticles const& /*mcparticles*/)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& mccollision : mccollisions) {
      auto collision_per_mccoll = collisions.sliceBy(rec_perMcCollision, mccollision.globalIndex());
      int nrec_per_mc = collision_per_mccoll.size();
      fRegistry.fill(HIST("Event/hNrecPerMCCollision"), nrec_per_mc);
    }

    for (auto& collision : collisions) {
      if ((pairtype == kPHOSPHOS || pairtype == kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue; // I don't know why this is necessary in simulation.
      }

      if (eventcuts.onlyKeepWeightedEvents && std::fabs(collision.weight() - 1.) < 1E-10) {
        continue;
      }

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto mccollision = collision.template emmcevent_as<TMCCollisions>();
      auto binned_data_pi0_gen = mccollision.generatedPi0();
      auto binned_data_eta_gen = mccollision.generatedEta();
      fillBinnedData<0>(binned_data_pi0_gen, collision.weight());
      fillBinnedData<1>(binned_data_eta_gen, collision.weight());
    } // end of collision loop
  }

  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Filter prefilter_pcm = ifnode(pcmcuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::v0photonkf::pfbderived == static_cast<uint16_t>(0), true);
  Filter prefilter_primaryelectron = ifnode(dileptoncuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::emprimaryelectron::pfbderived == static_cast<uint16_t>(0), true);

  void processAnalysis(FilteredMyCollisions const& collisions, MyMCCollisions const& mccollisions, aod::EMMCParticles const& mcparticles, Types const&... args)
  {
    if constexpr (pairtype == PairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      runTruePairing(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut, mccollisions, mcparticles);
      runGenInfo(collisions, mccollisions, mcparticles);
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runTruePairing(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut, mccollisions, mcparticles);
      runGenInfo(collisions, mccollisions, mcparticles);
    } else if constexpr (pairtype == PairType::kEMCEMC) {
      auto emcclusters = std::get<0>(std::tie(args...));
      runTruePairing(collisions, emcclusters, emcclusters, nullptr, nullptr, perCollision_emc, perCollision_emc, fEMCCut, fEMCCut, mccollisions, mcparticles);
      runGenInfo(collisions, mccollisions, mcparticles);
    }

    // else if constexpr (pairtype == PairType::kPHOSPHOS) {
    //   auto phosclusters = std::get<0>(std::tie(args...));
    //   runPairing(collisions, phosclusters, phosclusters, nullptr, nullptr, perCollision_phos, perCollision_phos, fPHOSCut, fPHOSCut, nullptr, nullptr);
    // }
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
  }
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processAnalysis, "process pair analysis", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processDummy, "Dummy function", true);
};
#endif // PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMAMC_H_
