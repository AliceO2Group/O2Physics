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
// This code runs loop over dalitz ee table for dalitz QC.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include <iterator>
#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "DCAFitter/DCAFitterN.h"
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
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;

using MyEMH = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrackWithCov>;

struct dielectronQC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgAnalysisType{"cfgAnalysisType", static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC), "kQC:0, kUPC:1, kFlowV2:2, kFlowV3:3, kFlowV4:4, kPolarization:5, kHFll:6"};
  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  // Configurable<bool> cfgDo_v2{"cfgDo_v2", false, "flag to analyze v2"};
  // Configurable<bool> cfgDo_v3{"cfgDo_v3", false, "flag to analyze v3"};
  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {VARIABLE_WIDTH, -M_PI / 2, -M_PI / 4, 0.0f, +M_PI / 4, +M_PI / 2}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  ConfigurableAxis ConfMeeBins{"ConfMeeBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00}, "mee bins for output histograms"};
  ConfigurableAxis ConfPteeBins{"ConfPteeBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00}, "pTee bins for output histograms"};
  ConfigurableAxis ConfDCAeeBins{"ConfDCAeeBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAee bins for output histograms"};

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

  DielectronCut fDielectonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dielectroncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.8, "max eta for single track"};
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

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  // o2::vertexing::DCAFitterN<2> fitter;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> cent_bin_edges;
  std::vector<float> zvtx_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;
  int nmod = -1; // this is for flow analysis

  float beamM1 = o2::constants::physics::MassProton; // mass of beam
  float beamM2 = o2::constants::physics::MassProton; // mass of beam
  float beamE1 = 0.f;                                // beam energy
  float beamE2 = 0.f;                                // beam energy
  float beamP1 = 0.f;                                // beam momentum
  float beamP2 = 0.f;                                // beam momentum

  void init(InitContext& /*context*/)
  {
    zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
    zvtx_bin_edges.erase(zvtx_bin_edges.begin());

    cent_bin_edges = std::vector<float>(ConfCentBins.value.begin(), ConfCentBins.value.end());
    cent_bin_edges.erase(cent_bin_edges.begin());

    ep_bin_edges = std::vector<float>(ConfEPBins.value.begin(), ConfEPBins.value.end());
    ep_bin_edges.erase(ep_bin_edges.begin());

    occ_bin_edges = std::vector<float>(ConfOccupancyBins.value.begin(), ConfOccupancyBins.value.end());
    occ_bin_edges.erase(occ_bin_edges.begin());

    emh_pos = new MyEMH(ndepth);
    emh_ele = new MyEMH(ndepth);

    DefineEMEventCut();
    DefineDileptonCut();
    addhistograms();

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // fitter.setPropagateToPCA(true);
    // fitter.setMaxR(5.f);
    // fitter.setMinParamChange(1e-3);
    // fitter.setMinRelChi2Change(0.9);
    // fitter.setMaxDZIni(1e9);
    // fitter.setMaxChi2(1e9);
    // fitter.setUseAbsDCA(true);
    // fitter.setWeightedFinalPCA(false);
    // fitter.setMatCorrType(matCorr);
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
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      mRunNumber = collision.runNumber();
      // fitter.setBz(d_bz);
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

  ~dielectronQC()
  {
    delete emh_pos;
    emh_pos = 0x0;
    delete emh_ele;
    emh_ele = 0x0;

    map_mixed_eventId_to_qvector.clear();

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();
  }

  void addhistograms()
  {
    // pair info
    const AxisSpec axis_mass{ConfMeeBins, "m_{ee} (GeV/c^{2})"};
    const AxisSpec axis_pt{ConfPteeBins, "p_{T,ee} (GeV/c)"};
    const AxisSpec axis_dca{ConfDCAeeBins, "DCA_{ee}^{3D} (#sigma)"};
    std::string_view qvec_det_names[3] = {"FT0M", "FT0A", "FT0C"};

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC)) {
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kUPC)) {
      const AxisSpec axis_aco{10, 0, 1.f, "#alpha = 1 - #frac{|#varphi_{l^{+}} - #varphi_{l^{-}}|}{#pi}"};
      const AxisSpec axis_asym_pt{10, 0, 1.f, "A = #frac{|p_{T,l^{+}} - p_{T,l^{-}}|}{|p_{T,l^{+}} + p_{T,l^{-}}|}"};
      const AxisSpec axis_dphi_e_ee{18, 0, M_PI, "#Delta#varphi = #varphi_{e} - #varphi_{ee} (rad.)"};
      const AxisSpec axis_cos_theta_cs{10, 0.f, 1.f, "|cos(#theta_{CS})|"};
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_aco, axis_asym_pt, axis_dphi_e_ee, axis_cos_theta_cs}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2)) {
      nmod = 2;
      const AxisSpec axis_sp2{20, -5.f, 5.f, Form("u_{2}^{ee} #upoint Q_{2}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_sp2}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");

      const AxisSpec axis_sp2_l1{20, -5.f, 5.f, Form("u_{2,e1} #upoint Q_{2,1}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      const AxisSpec axis_sp2_l2{20, -5.f, 5.f, Form("u_{2,e2} #upoint Q_{2,2}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      const AxisSpec axis_cos_2dphi1{10, -1.f, 1.f, "cos(2 (#varphi_{e1} - #varphi_{ee}))"};
      const AxisSpec axis_cos_2dphi2{10, -1.f, 1.f, "cos(2 (#varphi_{e2} - #varphi_{ee}))"};
      fRegistry.add("Pair/mix/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_sp2_l1, axis_cos_2dphi1, axis_sp2_l2, axis_cos_2dphi2}, true);
      fRegistry.add("Pair/mix/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lspp/");
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lsmm/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
      nmod = 3;
      const AxisSpec axis_sp3{20, -5.f, 5.f, Form("u_{3}^{ee} #upoint Q_{3}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_sp3}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");

      const AxisSpec axis_sp3_l1{20, -5.f, 5.f, Form("u_{3,e1} #upoint Q_{3,1}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      const AxisSpec axis_sp3_l2{20, -5.f, 5.f, Form("u_{3,e2} #upoint Q_{3,2}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      const AxisSpec axis_cos_3dphi1{10, -1.f, 1.f, "cos(3 (#varphi_{e1} - #varphi_{ee}))"};
      const AxisSpec axis_cos_3dphi2{10, -1.f, 1.f, "cos(3 (#varphi_{e2} - #varphi_{ee}))"};
      fRegistry.add("Pair/mix/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_sp3_l1, axis_cos_3dphi1, axis_sp3_l2, axis_cos_3dphi2}, true);
      fRegistry.add("Pair/mix/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lspp/");
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lsmm/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV4)) {
      nmod = 4;
      const AxisSpec axis_sp4{20, -5.f, 5.f, Form("u_{4}^{ee} #upoint Q_{4}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_sp4}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");

      const AxisSpec axis_sp4_l1{20, -5.f, 5.f, Form("u_{4,e1} #upoint Q_{4,1}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      const AxisSpec axis_sp4_l2{20, -5.f, 5.f, Form("u_{4,e2} #upoint Q_{4,2}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
      const AxisSpec axis_cos_4dphi1{10, -1.f, 1.f, "cos(4 (#varphi_{e1} - #varphi_{ee}))"};
      const AxisSpec axis_cos_4dphi2{10, -1.f, 1.f, "cos(4 (#varphi_{e2} - #varphi_{ee}))"};
      fRegistry.add("Pair/mix/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_sp4_l1, axis_cos_4dphi1, axis_sp4_l2, axis_cos_4dphi2}, true);
      fRegistry.add("Pair/mix/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lspp/");
      fRegistry.addClone("Pair/mix/uls/", "Pair/mix/lsmm/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kPolarization)) {
      const AxisSpec axis_cos_theta_cs{10, 0.f, 1.f, "|cos(#theta_{CS})|"};
      const AxisSpec axis_phi_cs{18, 0.f, M_PI, "|#varphi_{CS}| (rad.)"};
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_cos_theta_cs, axis_phi_cs}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kHFll)) {
      const AxisSpec axis_dphi_ee{18, 0, M_PI, "#Delta#varphi = #varphi_{e1} - #varphi_{e2} (rad.)"};
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_dphi_ee}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    } else { // same as kQC to avoid seg. fault
      fRegistry.add("Pair/same/uls/hs", "dielectron", kTHnSparseD, {axis_mass, axis_pt, axis_dca}, true);
      fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2D, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
      fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
      fRegistry.addClone("Pair/same/", "Pair/mix/");
    }

    // for track info
    fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/hDCA3DSigma", "DCA 3D;DCA_{3D} (#sigma);", kTH1F, {{100, 0.0f, 10.0f}}, false);
    fRegistry.add("Track/hDCAxy_Pt", "DCA_{xy} vs. pT;p_{T} (GeV/c);DCA_{xy} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hDCAz_Pt", "DCA_{z} vs. pT;p_{T} (GeV/c);DCA_{z} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFbeta", "TOF beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    fRegistry.add("Track/h1overTOFbeta", "TOF beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {1000, 0.8, 1.8}}, false);
    fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH1F, {{32, 0, 16}}, false);

    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms(&fRegistry, nmod > 0);
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

  void DefineDileptonCut()
  {
    fDielectonCut = DielectronCut("fDielectonCut", "fDielectonCut");

    // for pair
    fDielectonCut.SetMeeRange(dielectroncuts.cfg_min_mass, dielectroncuts.cfg_max_mass);
    fDielectonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dielectroncuts.cfg_phiv_intercept) / dielectroncuts.cfg_phiv_slope; });
    fDielectonCut.SetPairDCARange(dielectroncuts.cfg_min_pair_dca3d, dielectroncuts.cfg_max_pair_dca3d); // in sigma
    fDielectonCut.ApplyPhiV(dielectroncuts.cfg_apply_phiv);
    fDielectonCut.ApplyPrefilter(dielectroncuts.cfg_apply_pf);
    fDielectonCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectonCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);

    // for track
    fDielectonCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, 1e+10f);
    fDielectonCut.SetTrackEtaRange(-dielectroncuts.cfg_max_eta_track, +dielectroncuts.cfg_max_eta_track);
    fDielectonCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectonCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectonCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectonCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectonCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectonCut.SetMeanClusterSizeITSob(0, 16);
    fDielectonCut.SetMaxDcaXY(dielectroncuts.cfg_max_dcaxy);
    fDielectonCut.SetMaxDcaZ(dielectroncuts.cfg_max_dcaz);

    // for eID
    fDielectonCut.SetPIDScheme(dielectroncuts.cfg_pid_scheme);
    fDielectonCut.SetTPCNsigmaElRange(dielectroncuts.cfg_min_TPCNsigmaEl, dielectroncuts.cfg_max_TPCNsigmaEl);
    fDielectonCut.SetTPCNsigmaMuRange(dielectroncuts.cfg_min_TPCNsigmaMu, dielectroncuts.cfg_max_TPCNsigmaMu);
    fDielectonCut.SetTPCNsigmaPiRange(dielectroncuts.cfg_min_TPCNsigmaPi, dielectroncuts.cfg_max_TPCNsigmaPi);
    fDielectonCut.SetTPCNsigmaKaRange(dielectroncuts.cfg_min_TPCNsigmaKa, dielectroncuts.cfg_max_TPCNsigmaKa);
    fDielectonCut.SetTPCNsigmaPrRange(dielectroncuts.cfg_min_TPCNsigmaPr, dielectroncuts.cfg_max_TPCNsigmaPr);
    fDielectonCut.SetTOFNsigmaElRange(dielectroncuts.cfg_min_TOFNsigmaEl, dielectroncuts.cfg_max_TOFNsigmaEl);

    if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      o2::ml::OnnxModel* eid_bdt = new o2::ml::OnnxModel();
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

      fDielectonCut.SetPIDModel(eid_bdt);
    } // end of PID ML
  }

  template <int ev_id, typename TCollision, typename TTrack1, typename TTrack2>
  // bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, std::array<float ,2> q2bpos_mix, etc is better)
  bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, const std::array<float, 2> qvector_mix = {999.f, 999.f})
  {
    if constexpr (ev_id == 1) {
      if (t1.has_ambiguousElectrons() && t2.has_ambiguousElectrons()) {
        for (auto& possible_id1 : t1.ambiguousElectronsIds()) {
          for (auto& possible_id2 : t2.ambiguousElectronsIds()) {
            if (possible_id1 == possible_id2) {
              // LOGF(info, "event id = %d: same track is found. t1.trackId() = %d, t1.collisionId() = %d, t1.pt() = %f, t1.eta() = %f, t1.phi() = %f, t2.trackId() = %d, t2.collisionId() = %d, t2.pt() = %f, t2.eta() = %f, t2.phi() = %f", ev_id, t1.trackId(), t1.collisionId(), t1.pt(), t1.eta(), t1.phi(), t2.trackId(), t2.collisionId(), t2.pt(), t2.eta(), t2.phi());
              return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
            }
          }
        }
      }
    }

    if constexpr (ev_id == 0) {
      if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
        if (!fDielectonCut.IsSelectedTrack<true>(t1, collision) || !fDielectonCut.IsSelectedTrack<true>(t2, collision)) {
          return false;
        }
      } else { // cut-based
        if (!fDielectonCut.IsSelectedTrack(t1) || !fDielectonCut.IsSelectedTrack(t2)) {
          return false;
        }
      }
    }

    if (!fDielectonCut.IsSelectedPair(t1, t2, d_bz)) {
      return false;
    }

    // float pca = 999.f, lxy = 999.f; // in unit of cm
    // o2::aod::pwgem::dilepton::utils::pairutil::isSVFound(fitter, collision, t1, t2, pca, lxy);

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    if (abs(v12.Rapidity()) > maxY) {
      return false;
    }

    float dca_t1_3d = dca3DinSigma(t1);
    float dca_t2_3d = dca3DinSigma(t2);
    float dca_ee_3d = std::sqrt((dca_t1_3d * dca_t1_3d + dca_t2_3d * dca_t2_3d) / 2.);
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);

    //    if (t1.sign() * t2.sign() < 0) { // ULS
    //      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d);
    //      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
    //    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
    //      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d);
    //      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
    //    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
    //      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d);
    //      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
    //    }

    if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kQC)) {
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kUPC)) {
      float dphi = v1.Phi() - v2.Phi();
      o2::math_utils::bringToPMPi(dphi);
      float aco = 1.f - abs(dphi) / M_PI;
      float asym = abs(v1.Pt() - v2.Pt()) / (v1.Pt() + v2.Pt());
      float dphi_e_ee = v1.Phi() - v12.Phi();
      o2::math_utils::bringToPMPi(dphi_e_ee);

      float cos_thetaCS = 999, phiCS = 999.f;
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS<false>(t1, t2, o2::constants::physics::MassElectron, o2::constants::physics::MassElectron, beamE1, beamE2, beamP1, beamP2, cos_thetaCS, phiCS);

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d, aco, asym, abs(dphi_e_ee), abs(cos_thetaCS));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d, aco, asym, abs(dphi_e_ee), abs(cos_thetaCS));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d, aco, asym, abs(dphi_e_ee), abs(cos_thetaCS));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
      }

    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV2) || cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kFlowV3)) {
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
      std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};
      // const std::array<float, 2> q2vector[5] = {q2ft0m, q2ft0a, q2ft0c, q2bpos, q2bneg};
      std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
      std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
      std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
      std::array<float, 2> q3bpos = {collision.q3xbpos(), collision.q3ybpos()};
      std::array<float, 2> q3bneg = {collision.q3xbneg(), collision.q3ybneg()};
      // const std::array<float, 2> q3vector[5] = {q3ft0m, q3ft0a, q3ft0c, q3bpos, q3bneg};

      std::vector<std::vector<std::array<float, 2>>> qvectors = {
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 0th harmonics
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 1st harmonics
        {q2ft0m, q2ft0a, q2ft0c, q2bpos, q2bneg},                                         // 2nd harmonics
        {q3ft0m, q3ft0a, q3ft0c, q3bpos, q3bneg},                                         // 3rd harmonics
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 4th harmonics
      };

      float sp = 0.0;                 // for same event
      float sp1 = 999.f, sp2 = 999.f; // for mixed event
      float cos_dphi1 = 999.f, cos_dphi2 = 999.f;
      if constexpr (ev_id == 0) {
        sp = RecoDecay::dotProd(std::array<float, 2>{static_cast<float>(std::cos(nmod * v12.Phi())), static_cast<float>(std::sin(nmod * v12.Phi()))}, qvectors[nmod][cfgQvecEstimator]);

        if (t1.sign() * t2.sign() < 0) { // ULS
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d, sp);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
        } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d, sp);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
        } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d, sp);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
        }

      } else if constexpr (ev_id == 1) {
        sp1 = RecoDecay::dotProd(std::array<float, 2>{static_cast<float>(std::cos(nmod * v1.Phi())), static_cast<float>(std::sin(nmod * v1.Phi()))}, qvectors[nmod][cfgQvecEstimator]);
        sp2 = RecoDecay::dotProd(std::array<float, 2>{static_cast<float>(std::cos(nmod * v2.Phi())), static_cast<float>(std::sin(nmod * v2.Phi()))}, qvector_mix);
        cos_dphi1 = std::cos(nmod * (v1.Phi() - v12.Phi()));
        cos_dphi2 = std::cos(nmod * (v2.Phi() - v12.Phi()));

        if (t1.sign() * t2.sign() < 0) { // ULS
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d, sp1, cos_dphi1, sp2, cos_dphi2);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
        } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d, sp1, cos_dphi1, sp2, cos_dphi2);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
        } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d, sp1, cos_dphi1, sp2, cos_dphi2);
          fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
        }
      }
    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kPolarization)) {
      float cos_thetaCS = 999, phiCS = 999.f;
      o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS<false>(t1, t2, o2::constants::physics::MassElectron, o2::constants::physics::MassElectron, beamE1, beamE2, beamP1, beamP2, cos_thetaCS, phiCS);
      // LOGF(info, "cos_thetaCS = %f, phiCS = %f", cos_thetaCS, phiCS);
      o2::math_utils::bringToPMPi(phiCS);

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d, abs(cos_thetaCS), abs(phiCS));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d, abs(cos_thetaCS), abs(phiCS));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d, abs(cos_thetaCS), abs(phiCS));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
      }

    } else if (cfgAnalysisType == static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonAnalysisType::kHFll)) {
      float dphi = v1.Phi() - v2.Phi();
      o2::math_utils::bringToPMPi(dphi);

      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d, abs(dphi));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d, abs(dphi));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d, abs(dphi));
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
      }

    } else {                           // same as kQC to avoid seg. fault
      if (t1.sign() * t2.sign() < 0) { // ULS
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_ee_3d);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_ee_3d);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hMvsPhiV"), phiv, v12.M());
      } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_ee_3d);
        fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hMvsPhiV"), phiv, v12.M());
      }
    }

    // store tracks for event mixing without double counting
    if constexpr (ev_id == 0) {
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());
      std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, t1.globalIndex());
      std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, t2.globalIndex());

      std::vector<int> possibleIds1;
      std::vector<int> possibleIds2;
      std::copy(t1.ambiguousElectronsIds().begin(), t1.ambiguousElectronsIds().end(), std::back_inserter(possibleIds1));
      std::copy(t2.ambiguousElectronsIds().begin(), t2.ambiguousElectronsIds().end(), std::back_inserter(possibleIds2));

      if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id1) == used_trackIds.end()) {
        used_trackIds.emplace_back(pair_tmp_id1);
        fillTrackInfo(t1);
        if (cfgDoMix) {
          if (t1.sign() > 0) {
            emh_pos->AddTrackToEventPool(key_df_collision, EMTrackWithCov(t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1,
                                                                          t1.x(), t1.y(), t1.z(), t1.alpha(), t1.snp(), t1.tgl(), t1.cYY(), t1.cZY(), t1.cZZ(),
                                                                          t1.cSnpY(), t1.cSnpZ(), t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(), t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()));
          } else {
            emh_ele->AddTrackToEventPool(key_df_collision, EMTrackWithCov(t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1,
                                                                          t1.x(), t1.y(), t1.z(), t1.alpha(), t1.snp(), t1.tgl(), t1.cYY(), t1.cZY(), t1.cZZ(),
                                                                          t1.cSnpY(), t1.cSnpZ(), t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(), t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()));
          }
        }
      }
      if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id2) == used_trackIds.end()) {
        used_trackIds.emplace_back(pair_tmp_id2);
        fillTrackInfo(t2);
        if (cfgDoMix) {
          if (t2.sign() > 0) {
            emh_pos->AddTrackToEventPool(key_df_collision, EMTrackWithCov(t2.globalIndex(), collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron, t2.sign(), t2.dcaXY(), t2.dcaZ(), possibleIds1,
                                                                          t2.x(), t2.y(), t2.z(), t2.alpha(), t2.snp(), t2.tgl(), t2.cYY(), t2.cZY(), t2.cZZ(),
                                                                          t2.cSnpY(), t2.cSnpZ(), t2.cSnpSnp(), t2.cTglY(), t2.cTglZ(), t2.cTglSnp(), t2.cTglTgl(), t2.c1PtY(), t2.c1PtZ(), t2.c1PtSnp(), t2.c1PtTgl(), t2.c1Pt21Pt2()));
          } else {
            emh_ele->AddTrackToEventPool(key_df_collision, EMTrackWithCov(t2.globalIndex(), collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron, t2.sign(), t2.dcaXY(), t2.dcaZ(), possibleIds1,
                                                                          t2.x(), t2.y(), t2.z(), t2.alpha(), t2.snp(), t2.tgl(), t2.cYY(), t2.cZY(), t2.cZZ(),
                                                                          t2.cSnpY(), t2.cSnpZ(), t2.cSnpSnp(), t2.cTglY(), t2.cTglZ(), t2.cTglSnp(), t2.cTglTgl(), t2.c1PtY(), t2.c1PtZ(), t2.c1PtSnp(), t2.c1PtTgl(), t2.c1Pt21Pt2()));
          }
        }
      }
    }
    return true;
  }

  template <typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    float dca_3d = dca3DinSigma(track);
    fRegistry.fill(HIST("Track/hPt"), track.pt());
    fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/hEtaPhi"), track.phi(), track.eta());
    fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
    fRegistry.fill(HIST("Track/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
    fRegistry.fill(HIST("Track/hDCA3DSigma"), dca_3d);
    fRegistry.fill(HIST("Track/hDCAxy_Pt"), track.pt(), track.dcaXY());
    fRegistry.fill(HIST("Track/hDCAz_Pt"), track.pt(), track.dcaZ());
    fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
    fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
    fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
    fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
    fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
    fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
    fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
    fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
    fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
    fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
    fRegistry.fill(HIST("Track/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
    fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
    fRegistry.fill(HIST("Track/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
    fRegistry.fill(HIST("Track/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
    fRegistry.fill(HIST("Track/hTOFbeta"), track.tpcInnerParam(), track.beta());
    fRegistry.fill(HIST("Track/h1overTOFbeta"), track.tpcInnerParam(), 1. / track.beta());
    fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
    fRegistry.fill(HIST("Track/hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
    fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
    fRegistry.fill(HIST("Track/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
    fRegistry.fill(HIST("Track/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  SliceCache cache;
  Preslice<MyTracks> perCollision_track = aod::emprimaryelectron::emeventId;
  Filter trackFilter = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dielectroncuts.cfg_max_eta_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter = (dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < dielectroncuts.cfg_min_TPCNsigmaPi || dielectroncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi) && ((0.96f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.04f) || o2::aod::pidtofbeta::beta < 0.f);
  Filter ttcaFilter = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  using FilteredMyTracks = soa::Filtered<MyTracks>;
  using FilteredMyTrack = FilteredMyTracks::iterator;

  MyEMH* emh_pos = nullptr;
  MyEMH* emh_ele = nullptr;

  Partition<FilteredMyTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  std::map<std::pair<int, int>, std::array<float, 2>> map_mixed_eventId_to_qvector;

  std::vector<std::pair<int, int>> used_trackIds;
  int ndf = 0;
  void processQC(FilteredMyCollisions const& collisions, FilteredMyTracks const& /*tracks*/)
  {
    for (auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, nmod > 0);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, nmod > 0);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
      std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};
      // const std::array<float, 2> q2vector[5] = {q2ft0m, q2ft0a, q2ft0c, q2bpos, q2bneg};
      std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
      std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
      std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
      std::array<float, 2> q3bpos = {collision.q3xbpos(), collision.q3ybpos()};
      std::array<float, 2> q3bneg = {collision.q3xbneg(), collision.q3ybneg()};
      // const std::array<float, 2> q3vector[5] = {q3ft0m, q3ft0a, q3ft0c, q3bpos, q3bneg};

      std::vector<std::vector<std::array<float, 2>>> qvectors = {
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 0th harmonics
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 1st harmonics
        {q2ft0m, q2ft0a, q2ft0c, q2bpos, q2bneg},                                         // 2nd harmonics
        {q3ft0m, q3ft0a, q3ft0c, q3bpos, q3bneg},                                         // 3rd harmonics
        {{999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}, {999.f, 999.f}}, // 4th harmonics
      };

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      // LOGF(info, "collision.globalIndex() = %d , collision.posZ() = %f , collision.numContrib() = %d, centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", collision.globalIndex(), collision.posZ(), collision.numContrib(), centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      int nuls = 0, nlspp = 0, nlsmm = 0;
      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        bool is_pair_ok = fillPairInfo<0>(collision, pos, ele);
        if (is_pair_ok) {
          nuls++;
        }
      }
      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        bool is_pair_ok = fillPairInfo<0>(collision, pos1, pos2);
        if (is_pair_ok) {
          nlspp++;
        }
      }
      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        bool is_pair_ok = fillPairInfo<0>(collision, ele1, ele2);
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
      auto selected_negTracks_in_this_event = emh_ele->GetTracksPerCollision(key_df_collision);
      // LOGF(info, "N selected tracks in current event (%d, %d), zvtx = %f, centrality = %f , npos = %d , nele = %d, nuls = %d , nlspp = %d, nlsmm = %d", ndf, collision.globalIndex(), collision.posZ(), centralities[cfgCentEstimator], selected_posTracks_in_this_event.size(), selected_negTracks_in_this_event.size(), nuls, nlspp, nlsmm);

      auto collisionIds_in_mixing_pool = emh_pos->GetCollisionIdsFromEventPool(key_bin); // pos/ele does not matter.
      // LOGF(info, "collisionIds_in_mixing_pool.size() = %d", collisionIds_in_mixing_pool.size());

      for (auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
        int mix_dfId = mix_dfId_collisionId.first;
        int mix_collisionId = mix_dfId_collisionId.second;
        if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
          continue;
        }

        // float q2x_mixed = map_mixed_eventId_to_qvector[mix_dfId_collisionId][0];
        // float q2y_mixed = map_mixed_eventId_to_qvector[mix_dfId_collisionId][1];

        std::array<float, 2> qvector_mix = map_mixed_eventId_to_qvector[mix_dfId_collisionId];

        // LOGF(info, "mix_df = %d, mix_collisionId = %d, q2x_mixed = %f, q2y_mixed = %f", mix_dfId, mix_collisionId, q2x_mixed, q2y_mixed);

        auto posTracks_from_event_pool = emh_pos->GetTracksPerCollision(mix_dfId_collisionId);
        auto negTracks_from_event_pool = emh_ele->GetTracksPerCollision(mix_dfId_collisionId);
        // LOGF(info, "Do event mixing: current event (%d, %d) | event pool (%d, %d), npos = %d , nele = %d", ndf, collision.globalIndex(), mix_dfId, mix_collisionId, posTracks_from_event_pool.size(), negTracks_from_event_pool.size());

        for (auto& pos : selected_posTracks_in_this_event) { // ULS mix
          for (auto& ele : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos, ele, qvector_mix);
          }
        }

        for (auto& ele : selected_negTracks_in_this_event) { // ULS mix
          for (auto& pos : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, ele, pos, qvector_mix);
          }
        }

        for (auto& pos1 : selected_posTracks_in_this_event) { // LS++ mix
          for (auto& pos2 : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos1, pos2, qvector_mix);
          }
        }

        for (auto& ele1 : selected_negTracks_in_this_event) { // LS-- mix
          for (auto& ele2 : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, ele1, ele2, qvector_mix);
          }
        }
      } // end of loop over mixed event pool

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        // LOGF(info, "adding : ndf = %d, collision.globalIndex() = %d, q2x = %f, q2y  = %f", ndf, collision.globalIndex(), q2vector[cfgQvecEstimator][0], q2vector[cfgQvecEstimator][1]);
        map_mixed_eventId_to_qvector[key_df_collision] = qvectors[nmod][cfgQvecEstimator];
        emh_pos->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_ele->AddCollisionIdAtLast(key_bin, key_df_collision);
      }

    } // end of collision loop

    ndf++;
  } // end of process
  PROCESS_SWITCH(dielectronQC, processQC, "run dielectron QC", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(dielectronQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dielectronQC>(cfgc, TaskName{"dielectron-qc"})};
}
