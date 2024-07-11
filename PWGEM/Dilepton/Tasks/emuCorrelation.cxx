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

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EMFwdTrack.h"
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

using MyElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyElectron = MyElectrons::iterator;

using MyMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds>;
using MyMuon = MyMuons::iterator;

using MyEMH_el = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMTrackWithCov>;
using MyEMH_mu = o2::aod::pwgem::dilepton::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMFwdTrackWithCov>;

struct emuCorrelation {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 0, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 1000, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {VARIABLE_WIDTH, -M_PI / 2, -M_PI / 4, 0.0f, +M_PI / 4, +M_PI / 2}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  ConfigurableAxis ConfMemuBins{"ConfMemuBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0}, "memu bins for output histograms"};
  ConfigurableAxis ConfPtemuBins{"ConfPtemuBins", {VARIABLE_WIDTH, 0.00, 10.00}, "pTemu bins for output histograms"};
  ConfigurableAxis ConfDCAemuBins{"ConfDCAemuBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "single leg DCA bins for output histograms"};

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

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
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
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 80, "min ncrossed rows"};
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

  DimuonCut fDimuonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dimuoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pt"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pt"};
    Configurable<float> cfg_min_pair_dcaxy{"cfg_min_pair_dcaxy", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dcaxy{"cfg_max_pair_dcaxy", 1e+10, "max pair dca3d in sigma"};

    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+10, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 1e+10, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
  } dimuoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> cent_bin_edges;
  std::vector<float> zvtx_bin_edges;
  std::vector<float> ep_bin_edges;
  std::vector<float> occ_bin_edges;

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

    emh_el_pos = new MyEMH_el(ndepth);
    emh_mu_pos = new MyEMH_mu(ndepth);
    emh_el_neg = new MyEMH_el(ndepth);
    emh_mu_neg = new MyEMH_mu(ndepth);

    DefineEMEventCut();
    DefineDielectronCut();
    DefineDimuonCut();
    addhistograms();

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(5.f);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
    fitter.setMatCorrType(matCorr);
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
      fitter.setBz(d_bz);
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
    fitter.setBz(d_bz);
  }

  ~emuCorrelation()
  {
    delete emh_el_pos;
    emh_el_pos = 0x0;
    delete emh_mu_pos;
    emh_mu_pos = 0x0;

    delete emh_el_neg;
    emh_el_neg = 0x0;
    delete emh_mu_neg;
    emh_mu_neg = 0x0;

    used_electronIds.clear();
    used_electronIds.shrink_to_fit();
    used_muonIds.clear();
    used_muonIds.shrink_to_fit();
  }

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms(&fRegistry, false);

    // pair info
    const AxisSpec axis_mass{ConfMemuBins, "m_{e#mu} (GeV/c^{2})"};
    const AxisSpec axis_pt{ConfPtemuBins, "p_{T,e#mu} (GeV/c)"};
    const AxisSpec axis_rapidity{48, -4.0, +0.8, "y_{e#mu}"};
    const AxisSpec axis_dca_el{ConfDCAemuBins, "DCA_{e}^{3D} (#sigma)"};
    const AxisSpec axis_dca_mu{ConfDCAemuBins, "DCA_{#mu}^{XY} (#sigma)"};
    const AxisSpec axis_dphi{18, 0, M_PI, "#Delta#varphi_{e#mu} (rad.)"};

    // std::string_view qvec_det_names[3] = {"FT0M", "FT0A", "FT0C"};

    fRegistry.add("Pair/same/uls/hs", "e#mu", kTHnSparseD, {axis_mass, axis_pt, axis_rapidity, axis_dphi, axis_dca_el, axis_dca_mu}, true);
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
    fRegistry.addClone("Pair/same/", "Pair/mix/");

    // for electron info
    fRegistry.add("Track/Electron/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/Electron/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/Electron/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/Electron/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/Electron/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/Electron/hDCAxy_Pt", "DCA_{xy} vs. pT;p_{T} (GeV/c);DCA_{xy} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/Electron/hDCAz_Pt", "DCA_{z} vs. pT;p_{T} (GeV/c);DCA_{z} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/Electron/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/Electron/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/Electron/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/Electron/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/Electron/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/Electron/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/Electron/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTOFbeta", "TOF beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    fRegistry.add("Track/Electron/h1overTOFbeta", "TOF beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {1000, 0.8, 1.8}}, false);
    fRegistry.add("Track/Electron/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/Electron/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/Electron/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/Electron/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/Electron/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/Electron/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH1F, {{32, 0, 16}}, false);

    // for muon info
    fRegistry.add("Track/Muon/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/Muon/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/Muon/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {25, -4.5f, -2.0f}}, false);
    fRegistry.add("Track/Muon/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
    fRegistry.add("Track/Muon/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/Muon/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/Muon/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/Muon/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/Muon/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
    fRegistry.add("Track/Muon/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
    fRegistry.add("Track/Muon/hPDCA", "pDCA;p_{T} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/Muon/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/Muon/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/Muon/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/Muon/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
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

  void DefineDielectronCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for pair
    fDielectronCut.SetMeeRange(dielectroncuts.cfg_min_mass, dielectroncuts.cfg_max_mass);
    fDielectronCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dielectroncuts.cfg_phiv_intercept) / dielectroncuts.cfg_phiv_slope; });
    fDielectronCut.SetPairDCARange(dielectroncuts.cfg_min_pair_dca3d, dielectroncuts.cfg_max_pair_dca3d); // in sigma
    fDielectronCut.ApplyPhiV(dielectroncuts.cfg_apply_phiv);
    fDielectronCut.ApplyPrefilter(dielectroncuts.cfg_apply_pf);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, 1e+10f);
    fDielectronCut.SetTrackEtaRange(-dielectroncuts.cfg_max_eta_track, +dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITSob(0, 16);
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

      fDielectronCut.SetPIDModel(eid_bdt);
    } // end of PID ML
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // for pair
    fDimuonCut.SetMassRange(dimuoncuts.cfg_min_mass, dimuoncuts.cfg_max_mass);
    fDimuonCut.SetPairPtRange(dimuoncuts.cfg_min_pair_pt, dimuoncuts.cfg_max_pair_pt);
    fDimuonCut.SetPairDCAxyRange(dimuoncuts.cfg_min_pair_dcaxy, dimuoncuts.cfg_max_pair_dcaxy); // DCAxy in cm

    // for track
    fDimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, 1e10f);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 16);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
  }

  template <int ev_id, typename TCollision, typename TElectron, typename TMuon>
  bool fillPairInfo(TCollision const& collision, TElectron const& t1, TMuon const& t2)
  {
    // if constexpr (ev_id == 1) {
    //   if (t1.has_ambiguousElectrons() && t2.has_ambiguousMuons()) {
    //     for (auto& possible_id1 : t1.ambiguousElectronsIds()) {
    //       for (auto& possible_id2 : t2.ambiguousMuonsIds()) {
    //         if (possible_id1 == possible_id2) {
    //           // LOGF(info, "event id = %d: same track is found. t1.trackId() = %d, t1.collisionId() = %d, t1.pt() = %f, t1.eta() = %f, t1.phi() = %f, t2.trackId() = %d, t2.collisionId() = %d, t2.pt() = %f, t2.eta() = %f, t2.phi() = %f", ev_id, t1.trackId(), t1.collisionId(), t1.pt(), t1.eta(), t1.phi(), t2.trackId(), t2.collisionId(), t2.pt(), t2.eta(), t2.phi());
    //           return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions. This never happens in e-mu correlation.
    //         }
    //       }
    //     }
    //   }
    // }

    if constexpr (ev_id == 0) {
      if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
        if (!fDielectronCut.IsSelectedTrack<true>(t1, collision)) {
          return false;
        }
      } else { // cut-based
        if (!fDielectronCut.IsSelectedTrack(t1)) {
          return false;
        }
      }

      if (!fDimuonCut.IsSelectedTrack(t2)) {
        return false;
      }
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dphi = t1.phi() - t2.phi();
    if (-2 * M_PI < dphi && dphi < -M_PI) {
      dphi += 2 * M_PI;
    } else if (-M_PI < dphi && dphi < 0) {
      dphi += M_PI;
    } else if (dphi > M_PI) {
      dphi -= M_PI;
    }

    float dca_3d_el = dca3DinSigma(t1);
    float dca_xy_mu = fwdDcaXYinSigma(t2);

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
    }

    // store tracks for event mixing without double counting
    if constexpr (ev_id == 0) {
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());
      std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, t1.globalIndex());
      std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, t2.globalIndex());

      std::vector<int> possibleIds1;
      std::vector<int> possibleIds2;
      std::copy(t1.ambiguousElectronsIds().begin(), t1.ambiguousElectronsIds().end(), std::back_inserter(possibleIds1));
      std::copy(t2.ambiguousMuonsIds().begin(), t2.ambiguousMuonsIds().end(), std::back_inserter(possibleIds2));

      if (std::find(used_electronIds.begin(), used_electronIds.end(), pair_tmp_id1) == used_electronIds.end()) {
        used_electronIds.emplace_back(pair_tmp_id1);
        fillElectronInfo(t1);
        if (cfgDoMix) {
          if (t1.sign() > 0) {
            emh_el_pos->AddTrackToEventPool(key_df_collision, EMTrackWithCov(t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1,
                                                                             t1.x(), t1.y(), t1.z(), t1.alpha(), t1.snp(), t1.tgl(), t1.cYY(), t1.cZY(), t1.cZZ(),
                                                                             t1.cSnpY(), t1.cSnpZ(), t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(), t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()));
          } else {
            emh_el_neg->AddTrackToEventPool(key_df_collision, EMTrackWithCov(t1.globalIndex(), collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron, t1.sign(), t1.dcaXY(), t1.dcaZ(), possibleIds1,
                                                                             t1.x(), t1.y(), t1.z(), t1.alpha(), t1.snp(), t1.tgl(), t1.cYY(), t1.cZY(), t1.cZZ(),
                                                                             t1.cSnpY(), t1.cSnpZ(), t1.cSnpSnp(), t1.cTglY(), t1.cTglZ(), t1.cTglSnp(), t1.cTglTgl(), t1.c1PtY(), t1.c1PtZ(), t1.c1PtSnp(), t1.c1PtTgl(), t1.c1Pt21Pt2()));
          }
        }
      }
      if (std::find(used_muonIds.begin(), used_muonIds.end(), pair_tmp_id2) == used_muonIds.end()) {
        used_muonIds.emplace_back(pair_tmp_id2);
        fillMuonInfo(t2);
        if (cfgDoMix) {
          if (t2.sign() > 0) {
            emh_mu_pos->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                                t2.x(), t2.y(), t2.z(), t2.tgl(), t2.cXX(), t2.cXY(), t2.cYY(),
                                                                                t2.cPhiX(), t2.cPhiY(), t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(), t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2(), t2.chi2()));
          } else {
            emh_mu_neg->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                                t2.x(), t2.y(), t2.z(), t2.tgl(), t2.cXX(), t2.cXY(), t2.cYY(),
                                                                                t2.cPhiX(), t2.cPhiY(), t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(), t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2(), t2.chi2()));
          }
        }
      }
    }
    return true;
  }

  template <typename TTrack>
  void fillElectronInfo(TTrack const& track)
  {
    fRegistry.fill(HIST("Track/Electron/hPt"), track.pt());
    fRegistry.fill(HIST("Track/Electron/hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/Electron/hEtaPhi"), track.phi(), track.eta());
    fRegistry.fill(HIST("Track/Electron/hDCAxyz"), track.dcaXY(), track.dcaZ());
    fRegistry.fill(HIST("Track/Electron/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
    fRegistry.fill(HIST("Track/Electron/hDCAxy_Pt"), track.pt(), track.dcaXY());
    fRegistry.fill(HIST("Track/Electron/hDCAz_Pt"), track.pt(), track.dcaZ());
    fRegistry.fill(HIST("Track/Electron/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/Electron/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
    fRegistry.fill(HIST("Track/Electron/hNclsITS"), track.itsNCls());
    fRegistry.fill(HIST("Track/Electron/hNclsTPC"), track.tpcNClsFound());
    fRegistry.fill(HIST("Track/Electron/hNcrTPC"), track.tpcNClsCrossedRows());
    fRegistry.fill(HIST("Track/Electron/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("Track/Electron/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("Track/Electron/hChi2TPC"), track.tpcChi2NCl());
    fRegistry.fill(HIST("Track/Electron/hChi2ITS"), track.itsChi2NCl());
    fRegistry.fill(HIST("Track/Electron/hITSClusterMap"), track.itsClusterMap());
    fRegistry.fill(HIST("Track/Electron/hMeanClusterSizeITS"), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
    fRegistry.fill(HIST("Track/Electron/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
    fRegistry.fill(HIST("Track/Electron/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
    fRegistry.fill(HIST("Track/Electron/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
    fRegistry.fill(HIST("Track/Electron/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
    fRegistry.fill(HIST("Track/Electron/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
    fRegistry.fill(HIST("Track/Electron/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
    fRegistry.fill(HIST("Track/Electron/hTOFbeta"), track.tpcInnerParam(), track.beta());
    fRegistry.fill(HIST("Track/Electron/h1overTOFbeta"), track.tpcInnerParam(), 1. / track.beta());
    fRegistry.fill(HIST("Track/Electron/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
    fRegistry.fill(HIST("Track/Electron/hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
    fRegistry.fill(HIST("Track/Electron/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
    fRegistry.fill(HIST("Track/Electron/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
    fRegistry.fill(HIST("Track/Electron/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
  }

  template <typename TTrack>
  void fillMuonInfo(TTrack const& track)
  {
    fRegistry.fill(HIST("Track/Muon/hPt"), track.pt());
    fRegistry.fill(HIST("Track/Muon/hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/Muon/hEtaPhi"), track.phi(), track.eta());
    fRegistry.fill(HIST("Track/Muon/hTrackType"), track.trackType());
    fRegistry.fill(HIST("Track/Muon/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
    fRegistry.fill(HIST("Track/Muon/hDCAxySigma"), track.fwdDcaX() / sqrt(track.cXX()), track.fwdDcaY() / sqrt(track.cYY()));
    fRegistry.fill(HIST("Track/Muon/hDCAxRes_Pt"), track.pt(), sqrt(track.cXX()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/Muon/hDCAyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/Muon/hNclsMFT"), track.nClustersMFT());
    fRegistry.fill(HIST("Track/Muon/hNclsMCH"), track.nClusters());
    fRegistry.fill(HIST("Track/Muon/hPDCA"), track.pt(), track.pDca());
    fRegistry.fill(HIST("Track/Muon/hChi2"), track.chi2());
    fRegistry.fill(HIST("Track/Muon/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
    fRegistry.fill(HIST("Track/Muon/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
    fRegistry.fill(HIST("Track/Muon/hMFTClusterMap"), track.mftClusterMap());
  }
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  SliceCache cache;
  Preslice<MyElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Filter trackFilter_el = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dielectroncuts.cfg_max_eta_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter = (dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < dielectroncuts.cfg_min_TPCNsigmaPi || dielectroncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi) && ((0.96f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.04f) || o2::aod::pidtofbeta::beta < 0.f);
  Filter ttcaFilter_el = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  using FilteredMyElectrons = soa::Filtered<MyElectrons>;
  using FilteredMyElectron = FilteredMyElectrons::iterator;

  Preslice<MyMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_mu = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track;
  Filter ttcaFilter_mu = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);
  using FilteredMyMuons = soa::Filtered<MyMuons>;

  MyEMH_el* emh_el_pos = nullptr;
  MyEMH_mu* emh_mu_pos = nullptr;
  MyEMH_el* emh_el_neg = nullptr;
  MyEMH_mu* emh_mu_neg = nullptr;

  Partition<FilteredMyElectrons> posTracks_el = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyElectrons> negTracks_el = o2::aod::emprimaryelectron::sign < int8_t(0);
  Partition<FilteredMyMuons> posTracks_mu = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<FilteredMyMuons> negTracks_mu = o2::aod::emprimarymuon::sign < int8_t(0);

  std::vector<std::pair<int, int>> used_electronIds;
  std::vector<std::pair<int, int>> used_muonIds;
  int ndf = 0;
  void processQC(FilteredMyCollisions const& collisions, FilteredMyElectrons const&, FilteredMyMuons const&)
  {
    for (auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, false);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, false);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto posTracks_el_per_coll = posTracks_el->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_el_per_coll = negTracks_el->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto posTracks_mu_per_coll = posTracks_mu->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      auto negTracks_mu_per_coll = negTracks_mu->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);

      // LOGF(info, "posTracks_el_per_coll.size() = %d, negTracks_el_per_coll.size() = %d, posTracks_mu_per_coll.size() = %d, negTracks_mu_per_coll.size() = %d", posTracks_el_per_coll.size(), negTracks_el_per_coll.size(), posTracks_mu_per_coll.size(), negTracks_mu_per_coll.size() );

      int nuls = 0, nlspp = 0, nlsmm = 0;
      for (auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_el_per_coll, negTracks_mu_per_coll))) { // ULS+-
        bool is_pair_ok = fillPairInfo<0>(collision, pos, neg);
        if (is_pair_ok) {
          nuls++;
        }
      }

      for (auto& [neg, pos] : combinations(CombinationsFullIndexPolicy(negTracks_el_per_coll, posTracks_mu_per_coll))) { // ULS-+
        bool is_pair_ok = fillPairInfo<0>(collision, neg, pos);
        if (is_pair_ok) {
          nuls++;
        }
      }

      for (auto& [pos1, pos2] : combinations(CombinationsFullIndexPolicy(posTracks_el_per_coll, posTracks_mu_per_coll))) { // LS++
        bool is_pair_ok = fillPairInfo<0>(collision, pos1, pos2);
        if (is_pair_ok) {
          nlspp++;
        }
      }
      for (auto& [neg1, neg2] : combinations(CombinationsFullIndexPolicy(negTracks_el_per_coll, negTracks_mu_per_coll))) { // LS--
        bool is_pair_ok = fillPairInfo<0>(collision, neg1, neg2);
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
      auto selected_posTracks_el_in_this_event = emh_el_pos->GetTracksPerCollision(key_df_collision);
      auto selected_negTracks_el_in_this_event = emh_el_neg->GetTracksPerCollision(key_df_collision);
      auto selected_posTracks_mu_in_this_event = emh_mu_pos->GetTracksPerCollision(key_df_collision);
      auto selected_negTracks_mu_in_this_event = emh_mu_neg->GetTracksPerCollision(key_df_collision);

      auto collisionIds_in_mixing_pool = emh_el_pos->GetCollisionIdsFromEventPool(key_bin); // pos/ele does not matter.
      // LOGF(info, "collisionIds_in_mixing_pool.size() = %d", collisionIds_in_mixing_pool.size());

      for (auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
        int mix_dfId = mix_dfId_collisionId.first;
        int mix_collisionId = mix_dfId_collisionId.second;
        if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
          continue;
        }

        // LOGF(info, "mix_df = %d, mix_collisionId = %d, q2x_mixed = %f, q2y_mixed = %f", mix_dfId, mix_collisionId, q2x_mixed, q2y_mixed);

        auto posTracks_el_from_event_pool = emh_el_pos->GetTracksPerCollision(mix_dfId_collisionId);
        auto negTracks_el_from_event_pool = emh_el_neg->GetTracksPerCollision(mix_dfId_collisionId);
        auto posTracks_mu_from_event_pool = emh_mu_pos->GetTracksPerCollision(mix_dfId_collisionId);
        auto negTracks_mu_from_event_pool = emh_mu_neg->GetTracksPerCollision(mix_dfId_collisionId);
        // LOGF(info, "Do event mixing: current event (%d, %d) | event pool (%d, %d), npos = %d , nele = %d", ndf, collision.globalIndex(), mix_dfId, mix_collisionId, posTracks_from_event_pool.size(), negTracks_from_event_pool.size());

        for (auto& pos : selected_posTracks_el_in_this_event) { // ULS mix
          for (auto& neg : negTracks_mu_from_event_pool) {
            fillPairInfo<1>(collision, pos, neg);
          }
        }

        for (auto& neg : selected_negTracks_el_in_this_event) { // ULS mix
          for (auto& pos : posTracks_mu_from_event_pool) {
            fillPairInfo<1>(collision, neg, pos);
          }
        }

        for (auto& pos : selected_posTracks_mu_in_this_event) { // ULS mix
          for (auto& neg : negTracks_el_from_event_pool) {
            fillPairInfo<1>(collision, neg, pos);
          }
        }

        for (auto& neg : selected_negTracks_mu_in_this_event) { // ULS mix
          for (auto& pos : posTracks_el_from_event_pool) {
            fillPairInfo<1>(collision, pos, neg);
          }
        }

        for (auto& pos1 : selected_posTracks_el_in_this_event) { // LS++ mix
          for (auto& pos2 : posTracks_mu_from_event_pool) {
            fillPairInfo<1>(collision, pos1, pos2);
          }
        }

        for (auto& pos1 : selected_posTracks_mu_in_this_event) { // LS++ mix
          for (auto& pos2 : posTracks_el_from_event_pool) {
            fillPairInfo<1>(collision, pos2, pos1);
          }
        }

        for (auto& neg1 : selected_negTracks_el_in_this_event) { // LS-- mix
          for (auto& neg2 : negTracks_mu_from_event_pool) {
            fillPairInfo<1>(collision, neg1, neg2);
          }
        }

        for (auto& neg1 : selected_negTracks_mu_in_this_event) { // LS-- mix
          for (auto& neg2 : negTracks_el_from_event_pool) {
            fillPairInfo<1>(collision, neg2, neg1);
          }
        }

      } // end of loop over mixed event pool

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        // LOGF(info, "adding : ndf = %d, collision.globalIndex() = %d, q2x = %f, q2y  = %f", ndf, collision.globalIndex(), q2vector[cfgQvecEstimator][0], q2vector[cfgQvecEstimator][1]);
        emh_el_pos->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_mu_pos->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_el_neg->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_mu_neg->AddCollisionIdAtLast(key_bin, key_df_collision);
      }

    } // end of collision loop

    ndf++;
  } // end of process
  PROCESS_SWITCH(emuCorrelation, processQC, "run dielectron QC", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(emuCorrelation, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<emuCorrelation>(cfgc, TaskName{"emu-correlation"})};
}
