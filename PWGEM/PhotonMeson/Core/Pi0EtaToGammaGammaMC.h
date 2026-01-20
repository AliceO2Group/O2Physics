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
/// \file Pi0EtaToGammaGammaMC.h
/// \brief This code loops over photons and makes pairs for neutral mesons analyses with MC true info.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMAMC_H_
#define PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMAMC_H_

#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Utils/NMHistograms.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
// Dilepton headers
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
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

#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TH2.h>
#include <TString.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

template <o2::aod::pwgem::photonmeson::photonpair::PairType pairtype, o2::soa::is_table... Types>
struct Pi0EtaToGammaGammaMC {
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  o2::framework::Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  o2::framework::Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2"};
  o2::framework::Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  o2::framework::Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  o2::framework::Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  o2::framework::Configurable<float> maxY_rec{"maxY_rec", 0.9, "maximum rapidity for reconstructed particles"};
  o2::framework::Configurable<std::string> fd_k0s_to_pi0{"fd_k0s_pi0", "1.0", "feed down correction to pi0"};
  o2::framework::Configurable<bool> cfgRequireTrueAssociation{"cfgRequireTrueAssociation", false, "flag to require true mc collision association"};

  EMPhotonEventCut fEMEventCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "eventcut_group";
    o2::framework::Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    o2::framework::Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    o2::framework::Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
    o2::framework::Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    o2::framework::Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
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
    o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    o2::framework::Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    o2::framework::Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};

    o2::framework::Configurable<bool> cfg_apply_ml_cuts{"cfg_apply_ml", false, "flag to apply ML cut"};
    o2::framework::Configurable<bool> cfg_use_2d_binning{"cfg_use_2d_binning", true, "flag to use 2D binning (pT, cent)"};
    o2::framework::Configurable<bool> cfg_load_ml_models_from_ccdb{"cfg_load_ml_models_from_ccdb", true, "flag to load ML models from CCDB"};
    o2::framework::Configurable<int> cfg_timestamp_ccdb{"cfg_timestamp_ccdb", -1, "timestamp for CCDB"};
    o2::framework::Configurable<int> cfg_nclasses_ml{"cfg_nclasses_ml", static_cast<int>(o2::analysis::em_cuts_ml::NCutScores), "number of classes for ML"};
    o2::framework::Configurable<std::vector<int>> cfg_cut_dir_ml{"cfg_cut_dir_ml", std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}, "cut direction for ML"};
    o2::framework::Configurable<std::vector<std::string>> cfg_input_feature_names{"cfg_input_feature_names", std::vector<std::string>{"feature1", "feature2"}, "input feature names for ML models"};
    o2::framework::Configurable<std::vector<std::string>> cfg_model_paths_ccdb{"cfg_model_paths_ccdb", std::vector<std::string>{"path_ccdb/BDT_PCM/"}, "CCDB paths for ML models"};
    o2::framework::Configurable<std::vector<std::string>> cfg_onnx_file_names{"cfg_onnx_file_names", std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}, "ONNX file names for ML models"};
    o2::framework::Configurable<std::vector<double>> cfg_bins_pt_ml{"cfg_bins_pt_ml", std::vector<double>{o2::analysis::em_cuts_ml::vecBinsPt}, "pT bins for ML"};
    o2::framework::Configurable<std::vector<double>> cfg_bins_cent_ml{"cfg_bins_cent_ml", std::vector<double>{o2::analysis::em_cuts_ml::vecBinsCent}, "centrality bins for ML"};
    o2::framework::Configurable<o2::framework::LabeledArray<double>> cfg_cuts_pcm_ml{"cfg_cuts_pcm_ml", {o2::analysis::em_cuts_ml::Cuts[0], o2::analysis::em_cuts_ml::NBinsPt, o2::analysis::em_cuts_ml::NCutScores, o2::analysis::em_cuts_ml::labelsPt, o2::analysis::em_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
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
    o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    o2::framework::Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.05, "max dca XY for single track in cm"};
    o2::framework::Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.05, "max dca Z for single track in cm"};
    o2::framework::Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1.5, "max DCA 3D in sigma"};
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
    o2::framework::Configurable<bool> emcUseTM{"emcUseTM", true, "flag to use EMCal track matching cut or not."};
    o2::framework::Configurable<bool> emcUseSecondaryTM{"emcUseSecondaryTM", false, "flag to use EMCal secondary track matching cut or not. Required for PCM-EMC analyses."};
  } emccuts;

  o2::framework::Configurable<float> maxY_gen{"maxY_gen", 0.9, "maximum rapidity for generated particles"}; // for PCM and dielectron
  o2::framework::Configurable<float> maxRgen{"maxRgen", 90.f, "maximum radius for generated particles"};
  o2::framework::Configurable<float> margin_z_mc{"margin_z_mc", 7.0, "margin for z cut in cm for MC"};

  PHOSPhotonCut fPHOSCut;
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "phoscut_group";
    o2::framework::Configurable<float> cfg_min_Ecluster{"cfg_min_Ecluster", 0.3, "Minimum cluster energy for PHOS in GeV"};
  } phoscuts;

  TF1* f1fd_k0s_to_pi0;
  o2::framework::HistogramRegistry fRegistry{"output", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject, false, false};
  // static constexpr std::string_view event_types[2] = {"before/", "after/"};
  // static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};
  static constexpr std::string_view kParnames[2] = {"Pi0/", "Eta/"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  void init(o2::framework::InitContext&)
  {
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {
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

    fRegistry.add("Event/hNrecPerMCCollision", "Nrec per mc collision;N_{rec} collisions per MC collision", o2::framework::kTH1F, {{21, -0.5f, 20.5f}}, false);

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
    fV0PhotonCut.SetCutDirMl(pcmcuts.cfg_cut_dir_ml);
    fV0PhotonCut.SetMlModelPathsCCDB(pcmcuts.cfg_model_paths_ccdb);
    fV0PhotonCut.SetMlOnnxFileNames(pcmcuts.cfg_onnx_file_names);
    fV0PhotonCut.SetBinsPtMl(pcmcuts.cfg_bins_pt_ml);
    fV0PhotonCut.SetBinsCentMl(pcmcuts.cfg_bins_cent_ml);
    fV0PhotonCut.SetCutsPCMMl(pcmcuts.cfg_cuts_pcm_ml);
    fV0PhotonCut.SetNamesInputFeatures(pcmcuts.cfg_input_feature_names);

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

  o2::framework::SliceCache cache;
  o2::framework::PresliceOptional<o2::soa::Filtered<o2::soa::Join<o2::aod::V0PhotonsKF, o2::aod::V0KFEMEventIds, o2::aod::V0PhotonsKFPrefilterBitDerived>>> perCollision_pcm = o2::aod::v0photonkf::emeventId;
  o2::framework::PresliceOptional<o2::soa::Join<o2::aod::SkimEMCClusters, o2::aod::EMEMCClusterMCLabels, o2::aod::EMCEMEventIds>> perCollision_emc = o2::aod::emccluster::emeventId;
  o2::framework::PresliceOptional<o2::soa::Join<o2::aod::PHOSClusters, o2::aod::PHOSEMEventIds>> perCollision_phos = o2::aod::phoscluster::emeventId;
  o2::framework::PresliceOptional<o2::soa::Filtered<o2::soa::Join<o2::aod::EMPrimaryElectronsFromDalitz, o2::aod::EMPrimaryElectronEMEventIds, o2::aod::EMPrimaryElectronsPrefilterBitDerived, o2::aod::EMPrimaryElectronMCLabels>>> perCollision_electron = o2::aod::emprimaryelectron::emeventId;

  o2::framework::PresliceOptional<o2::aod::EmEmcMTracks> perEMCClusterMT = o2::aod::trackmatching::emEmcClusterId;
  o2::framework::PresliceOptional<o2::aod::EmEmcMSTracks> perEMCClusterMS = o2::aod::trackmatching::emEmcClusterId;

  o2::framework::Partition<o2::soa::Filtered<o2::soa::Join<o2::aod::EMPrimaryElectronsFromDalitz, o2::aod::EMPrimaryElectronEMEventIds, o2::aod::EMPrimaryElectronsPrefilterBitDerived, o2::aod::EMPrimaryElectronMCLabels>>> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt&& nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;
  o2::framework::Partition<o2::soa::Filtered<o2::soa::Join<o2::aod::EMPrimaryElectronsFromDalitz, o2::aod::EMPrimaryElectronEMEventIds, o2::aod::EMPrimaryElectronsPrefilterBitDerived, o2::aod::EMPrimaryElectronMCLabels>>> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;

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

    template <typename Self, typename Photon>
    static bool applyCut(Self& s, Photon const& g)
    {
      return s.fV0PhotonCut.template IsSelected<decltype(g), o2::soa::Join<o2::aod::V0Legs, o2::aod::V0LegMCLabels>>(g);
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
    template <typename Self, typename Cluster, typename TMatchedTracks = std::nullptr_t, typename TMatchedSecondaries = std::nullptr_t>
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

    template <typename Self, typename Cluster>
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
    template <typename Self, typename Track1, typename Track2>
    static bool applyCut(Self& s, Track1 const& track1, Track2 const& track2, float bz)
    {
      return s.fDileptonCut.IsSelected(track1, track2, bz);
    }
  };

  /// \brief function to run the photon pairing
  /// \tparam TDetectorTag1 tag for TPhotons1 type to select the proper cut function and arguments
  /// \tparam TDetectorTag2 tag for TPhotons2 type to select the proper cut function and arguments
  /// \tparam TCombinationPolicy StrictlyUpperIndex or FullIndex depending if we combine from the same detector or different detectors
  /// \tparam TCollisions collision table type
  /// \tparam TPhotons1 first photon table type
  /// \tparam TPhotons2 second photon table type
  /// \tparam TMCCollisions mc collision table type
  /// \tparam TMCParticles mc particles table type
  /// \tparam TLeg V0 leg table type used in TCut1 and/or TCut2 (for PCM/DalitzEE)
  /// \tparam TMatched matched global track table type for EMCal (for EMC)
  /// \tparam TSecondaries matched secondary track table type for EMCal (for EMC)
  /// \param collisions table of collisions (TCollisions)
  /// \param photons1 table of photons (TPhotons1)
  /// \param photons2 table of photons (TPhotons2)
  /// \param mccollisions table of mc collisions (TMCCollisions)
  /// \param mcparticles table of mc particles (TMCParticles)
  /// \param legs placeholder argument used only for template deduction (for PCM/DalitzEE)
  /// \param matchedtracks table of matched global tracks to EMCal clusters (for EMC)
  /// \param matchedsecondaries table of matched secondary tracks to EMCal clusters (for EMC)
  template <typename TDetectorTag1, typename TDetectorTag2, template <typename...> class TCombinationPolicy = o2::soa::CombinationsStrictlyUpperIndexPolicy, o2::soa::is_table TCollisions, o2::soa::is_table TPhotons1, o2::soa::is_table TPhotons2, o2::soa::is_table TMCCollisions, o2::soa::is_table TMCParticles, typename TLegs = std::nullptr_t, typename TMatchedTracks = std::nullptr_t, typename TMatchedSecondaries = std::nullptr_t>
  void runTruePairing(TCollisions const& collisions,
                      TPhotons1 const& photons1, TPhotons2 const& photons2,
                      TMCCollisions const& mccollisions, TMCParticles const& mcparticles,
                      TLegs const& /*legs*/ = nullptr, TMatchedTracks const& matchedTracks = nullptr, TMatchedSecondaries const& matchedSecondaries = nullptr)
  {
    for (auto& collision : collisions) {
      initCCDB(collision);
      fV0PhotonCut.SetCentrality(collision.centFT0M());
      fV0PhotonCut.SetD_Bz(d_bz);
      if ((pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS || pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }

      float weight = 1.f;
      if constexpr (std::is_same_v<std::decay_t<TCollisions>, o2::soa::Filtered<o2::soa::Join<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMMCEventLabels>, o2::aod::EMEventsWeight>>>) {
        weight = collision.weight();
      }

      if (eventcuts.onlyKeepWeightedEvents && std::fabs(weight - 1.0) < 1e-10) {
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

      int photonid1 = -1, photonid2 = -1, pi0id = -1, etaid = -1;
      if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM || pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS || pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) { // same kinds pairing

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

          if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM) { // check 2 legs
            auto pos1 = g1.template posTrack_as<TLegs>();
            auto ele1 = g1.template negTrack_as<TLegs>();
            auto pos2 = g2.template posTrack_as<TLegs>();
            auto ele2 = g2.template negTrack_as<TLegs>();

            auto pos1mc = pos1.template emmcparticle_as<TMCParticles>();
            auto ele1mc = ele1.template emmcparticle_as<TMCParticles>();
            auto pos2mc = pos2.template emmcparticle_as<TMCParticles>();
            auto ele2mc = ele2.template emmcparticle_as<TMCParticles>();

            photonid1 = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
            photonid2 = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles);
          } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) {
            auto cluster1mcparticle = mcparticles.iteratorAt(g1.emmcparticleId());
            auto cluster2mcparticle = mcparticles.iteratorAt(g2.emmcparticleId());

            photonid1 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(cluster1mcparticle, mcparticles, std::vector<int>{111, 221});
            photonid2 = o2::aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(cluster2mcparticle, mcparticles, std::vector<int>{111, 221});
          } else {
            photonid1 = -1;
            photonid2 = -1;
          }

          if (photonid1 < 0 || photonid2 < 0) {
            continue;
          }
          auto g1mc = mcparticles.iteratorAt(photonid1);
          auto g2mc = mcparticles.iteratorAt(photonid2);

          if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM) {
            if (!o2::aod::pwgem::photonmeson::utils::mcutil::IsConversionPointInAcceptance(g1mc, maxRgen, maxY_gen, margin_z_mc, mcparticles) || !o2::aod::pwgem::photonmeson::utils::mcutil::IsConversionPointInAcceptance(g2mc, maxRgen, maxY_gen, margin_z_mc, mcparticles)) {
              continue;
            }
          }

          pi0id = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
          etaid = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 221, mcparticles);

          if (g1mc.globalIndex() != g2mc.globalIndex() && pi0id < 0 && etaid < 0) { // for same gamma no pi0/eta will be found, but we still want to fill the FromSameGamma hist
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (std::fabs(v12.Rapidity()) > maxY_rec) {
            continue;
          }

          if (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) {
            float openingAngle = std::acos(v1.Vect().Dot(v2.Vect()) / (v1.P() * v2.P()));
            if (openingAngle < emccuts.minOpenAngle) {
              continue;
            }
          }

          if (g1mc.globalIndex() == g2mc.globalIndex()) {
            if (o2::aod::pwgem::dilepton::utils::mcutil::getMotherPDGCode(g1mc, mcparticles) == 111)
              fRegistry.fill(HIST("Pair/Pi0/hs_FromSameGamma"), v12.M(), v12.Pt(), weight);
            else if (o2::aod::pwgem::dilepton::utils::mcutil::getMotherPDGCode(g1mc, mcparticles) == 221)
              fRegistry.fill(HIST("Pair/Eta/hs_FromSameGamma"), v12.M(), v12.Pt(), weight);
            continue;
          }

          if (pi0id > 0) {
            auto pi0mc = mcparticles.iteratorAt(pi0id);
            if (cfgRequireTrueAssociation && (pi0mc.emmceventId() != collision.emmceventId())) {
              continue;
            }
            o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, pi0mc, mcparticles, mccollisions, f1fd_k0s_to_pi0, weight);
          } else if (etaid > 0) {
            auto etamc = mcparticles.iteratorAt(etaid);
            if (cfgRequireTrueAssociation && (etamc.emmceventId() != collision.emmceventId())) {
              continue;
            }
            o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, etamc, mcparticles, mccollisions, f1fd_k0s_to_pi0, weight);
          }
        } // end of pairing loop
      } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {

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
          auto pos1mc = pos1.template emmcparticle_as<TMCParticles>();
          auto ele1mc = ele1.template emmcparticle_as<TMCParticles>();
          photonid1 = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
          if (photonid1 < 0) {
            continue;
          }
          auto g1mc = mcparticles.iteratorAt(photonid1);
          if (!o2::aod::pwgem::photonmeson::utils::mcutil::IsConversionPointInAcceptance(g1mc, maxRgen, maxY_gen, margin_z_mc, mcparticles)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.f);

          for (const auto& [pos2, ele2] : o2::soa::combinations(TCombinationPolicy(positrons_per_collision, electrons_per_collision))) { // ULS
            if (pos2.trackId() == ele2.trackId()) {                                                                                      // this is protection against pairing identical 2 tracks.
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

            auto pos2mc = mcparticles.iteratorAt(pos2.emmcparticleId());
            auto ele2mc = mcparticles.iteratorAt(ele2.emmcparticleId());
            pi0id = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 111, mcparticles);
            etaid = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 221, mcparticles);
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
              o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, veeg, pi0mc, mcparticles, mccollisions, f1fd_k0s_to_pi0, weight);
            } else if (etaid > 0) {
              auto etamc = mcparticles.iteratorAt(etaid);
              if (cfgRequireTrueAssociation && (etamc.emmceventId() != collision.emmceventId())) {
                continue;
              }
              o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, veeg, etamc, mcparticles, mccollisions, f1fd_k0s_to_pi0, weight);
            }
          } // end of dielectron loop
        } // end of pcm loop
      } else { // PCM-EMC, PCM-PHOS.
        // TODO: implement proper functionality if we ever want to run this in Pb-Pb
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
          if (std::fabs(v12.Rapidity()) > maxY_rec) {
            continue;
          }
          // if (pi0id > 0) {
          //   auto pi0mc = mcparticles.iteratorAt(pi0id);
          //   o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, pi0mc, mcparticles, mccollisions, f1fd_k0s_to_pi0, weight);
          // } else if (etaid > 0) {
          //   auto etamc = mcparticles.iteratorAt(etaid);
          //   o2::aod::pwgem::photonmeson::utils::nmhistogram::fillTruePairInfo(&fRegistry, v12, etamc, mcparticles, mccollisions, f1fd_k0s_to_pi0, weight);
          // }
        } // end of pairing loop
      } // end of pairing in same event
    } // end of collision loop
  }

  template <int par_id, typename TBinnedData>
  void fillBinnedData(TBinnedData const& binned_data, const float weight = 1.f)
  {
    int xbin = 0, ybin = 0, zbin = 0;
    auto hPtY = fRegistry.get<TH2>(HIST("Generated/") + HIST(kParnames[par_id]) + HIST("hPtY")); // 2D
    auto hPt = fRegistry.get<TH1>(HIST("Generated/") + HIST(kParnames[par_id]) + HIST("hPt"));   // 1D

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

  o2::framework::PresliceUnsorted<o2::aod::EMMCParticles> perMcCollision = o2::aod::emmcparticle::emmceventId;
  o2::framework::PresliceUnsorted<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMMCEventLabels>> rec_perMcCollision = o2::aod::emmceventlabel::emmceventId;

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
      if ((pairtype == o2::aod::pwgem::photonmeson::photonpair::kPHOSPHOS || pairtype == o2::aod::pwgem::photonmeson::photonpair::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue; // I don't know why this is necessary in simulation.
      }

      float weight = 1.f;
      if constexpr (std::is_same_v<std::decay_t<TCollisions>, o2::soa::Filtered<o2::soa::Join<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMMCEventLabels>, o2::aod::EMEventsWeight>>>) {
        weight = collision.weight();
      }

      if (eventcuts.onlyKeepWeightedEvents && std::fabs(weight - 1.0) < 1e-10) {
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
      fillBinnedData<0>(binned_data_pi0_gen, weight);
      fillBinnedData<1>(binned_data_eta_gen, weight);
    } // end of collision loop
  }

  o2::framework::expressions::Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  o2::framework::expressions::Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  o2::framework::expressions::Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  // using FilteredMyCollisions = o2::soa::Filtered<o2::soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>>;

  o2::framework::expressions::Filter prefilter_pcm = ifnode(pcmcuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::v0photonkf::pfbderived == static_cast<uint16_t>(0), true);
  o2::framework::expressions::Filter prefilter_primaryelectron = ifnode(dileptoncuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::emprimaryelectron::pfbderived == static_cast<uint16_t>(0), true);

  void processAnalysis(o2::soa::Filtered<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMMCEventLabels>> const& collisions, o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts> const& mccollisions, o2::aod::EMMCParticles const& mcparticles, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM) {
      auto&& [v0photons, v0legs] = std::forward_as_tuple(args...);
      runTruePairing<PCMTag, PCMTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, v0photons, v0photons, mccollisions, mcparticles, v0legs);
      runGenInfo(collisions, mccollisions, mcparticles);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {
      auto&& [v0photons, v0legs, emprimaryelectrons] = std::forward_as_tuple(args...);
      runTruePairing<PCMTag, DalitzEETag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, emprimaryelectrons, mccollisions, mcparticles, v0legs);
      runGenInfo(collisions, mccollisions, mcparticles);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) {
      auto&& [emcclusters, emcMatchedTracks, emcMatchedSecondaries] = std::forward_as_tuple(args...);
      runTruePairing<EMCTag, EMCTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, emcclusters, emcclusters, mccollisions, mcparticles, nullptr, emcMatchedTracks, emcMatchedSecondaries);
      runGenInfo(collisions, mccollisions, mcparticles);
    }

    // else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS) {
    //   auto phosclusters = std::get<0>(std::tie(args...));
    //   runPairing(collisions, phosclusters, phosclusters, nullptr, nullptr, perCollision_phos, perCollision_phos, fPHOSCut, fPHOSCut, nullptr, nullptr);
    // }
    // else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMEMC) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto emcclusters = std::get<2>(std::tie(args...));
    //   auto emcmatchedtracks = std::get<3>(std::tie(args...));
    //   runPairing(collisions, v0photons, emcclusters, v0legs, nullptr, perCollision_pcm, perCollision_emc, fV0PhotonCut, fEMCCut, emcmatchedtracks, nullptr);
    // } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPHOS) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto phosclusters = std::get<2>(std::tie(args...));
    //   runPairing(collisions, v0photons, phosclusters, v0legs, nullptr, perCollision_pcm, perCollision_phos, fV0PhotonCut, fPHOSCut, nullptr, nullptr);
    // }
  }
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processAnalysis, "process pair analysis", true);

  // using FilteredMyCollisionsWithJJMC = o2::soa::Filtered<o2::soa::Join<o2::soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>, aod::EMEventsWeight>>;
  void processAnalysisJJMC(o2::soa::Filtered<o2::soa::Join<o2::soa::Join<o2::aod::EMEvents, o2::aod::EMEventsAlias, o2::aod::EMEventsMult, o2::aod::EMEventsCent, o2::aod::EMMCEventLabels>, o2::aod::EMEventsWeight>> const& collisions, o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts> const& mccollisions, o2::aod::EMMCParticles const& mcparticles, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPCM) {
      auto&& [v0photons, v0legs] = std::forward_as_tuple(args...);
      // auto v0photons = std::get<0>(std::tie(args...));
      // auto v0legs = std::get<1>(std::tie(args...));
      // runTruePairing(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut, mccollisions, mcparticles);
      runTruePairing<PCMTag, PCMTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, v0photons, v0photons, mccollisions, mcparticles, v0legs);
      runGenInfo(collisions, mccollisions, mcparticles);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMDalitzEE) {
      auto&& [v0photons, v0legs, emprimaryelectrons] = std::forward_as_tuple(args...);
      // auto v0photons = std::get<0>(std::tie(args...));
      // auto v0legs = std::get<1>(std::tie(args...));
      // auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      // runTruePairing(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut, mccollisions, mcparticles);
      runTruePairing<PCMTag, DalitzEETag, o2::soa::CombinationsFullIndexPolicy>(collisions, v0photons, emprimaryelectrons, mccollisions, mcparticles, v0legs);
      runGenInfo(collisions, mccollisions, mcparticles);
    } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kEMCEMC) {
      auto&& [emcclusters, emcMatchedTracks, emcMatchedSecondaries] = std::forward_as_tuple(args...);
      // auto emcclusters = std::get<0>(std::tie(args...));
      // runTruePairing(collisions, emcclusters, emcclusters, nullptr, nullptr, perCollision_emc, perCollision_emc, fEMCCut, fEMCCut, mccollisions, mcparticles);
      runTruePairing<EMCTag, EMCTag, o2::soa::CombinationsStrictlyUpperIndexPolicy>(collisions, emcclusters, emcclusters, mccollisions, mcparticles, nullptr, emcMatchedTracks, emcMatchedSecondaries);
      runGenInfo(collisions, mccollisions, mcparticles);
    }

    // else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPHOSPHOS) {
    //   auto phosclusters = std::get<0>(std::tie(args...));
    //   runPairing(collisions, phosclusters, phosclusters, nullptr, nullptr, perCollision_phos, perCollision_phos, fPHOSCut, fPHOSCut, nullptr, nullptr);
    // }
    // else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMEMC) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto emcclusters = std::get<2>(std::tie(args...));
    //   auto emcmatchedtracks = std::get<3>(std::tie(args...));
    //   runPairing(collisions, v0photons, emcclusters, v0legs, nullptr, perCollision_pcm, perCollision_emc, fV0PhotonCut, fEMCCut, emcmatchedtracks, nullptr);
    // } else if constexpr (pairtype == o2::aod::pwgem::photonmeson::photonpair::PairType::kPCMPHOS) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto phosclusters = std::get<2>(std::tie(args...));
    //   runPairing(collisions, v0photons, phosclusters, v0legs, nullptr, perCollision_pcm, perCollision_phos, fV0PhotonCut, fPHOSCut, nullptr, nullptr);
    // }
  }
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processAnalysisJJMC, "process pair analysis", false);

  void processDummy(o2::aod::EMEvents const&) {}
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processDummy, "Dummy function", false);
};
#endif // PWGEM_PHOTONMESON_CORE_PI0ETATOGAMMAGAMMAMC_H_
