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
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit, aod::EMPrimaryElectronMCLabels>;
using MyElectron = MyElectrons::iterator;

using MyMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMPrimaryMuonMCLabels>;
using MyMuon = MyMuons::iterator;

struct emuCorrelationMC {

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

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt_electron{"min_mcPt_electron", 0.05, "min. MC pT"};
    Configurable<float> max_mcPt_electron{"max_mcPt_electron", 1e+10, "max. MC pT"};
    Configurable<float> min_mcEta_electron{"min_mcEta_electron", -0.9, "max. MC eta"};
    Configurable<float> max_mcEta_electron{"max_mcEta_electron", +0.9, "max. MC eta"};

    Configurable<float> min_mcPt_muon{"min_mcPt_muon", 0.05, "min. MC pT"};
    Configurable<float> max_mcPt_muon{"max_mcPt_muon", 1e+10, "max. MC pT"};
    Configurable<float> min_mcEta_muon{"min_mcEta_muon", -4.0, "min. MC eta"};
    Configurable<float> max_mcEta_muon{"max_mcEta_muon", -2.5, "max. MC eta"};
  } mctrackcuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view lepton_source_types[6] = {"c2e/", "b2e/", "b2c2e/", "c2mu/", "b2mu/", "b2c2mu/"};

  void init(InitContext& /*context*/)
  {
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

  ~emuCorrelationMC() {}

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

    fRegistry.add("Pair/ccbar/c2e_c2mu/hadron_hadron/hs", "hs pair", kTHnSparseF, {axis_mass, axis_pt, axis_rapidity, axis_dphi, axis_dca_el, axis_dca_mu}, true);
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/hadron_hadron/", "Pair/ccbar/c2e_c2mu/meson_meson/");
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/hadron_hadron/", "Pair/ccbar/c2e_c2mu/baryon_baryon/");
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/hadron_hadron/", "Pair/ccbar/c2e_c2mu/meson_baryon/");
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/", "Pair/bbbar/b2e_b2mu/");
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/", "Pair/bbbar/b2c2e_b2c2mu/");
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/", "Pair/bbbar/b2c2e_b2mu_sameb/");
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/", "Pair/bbbar/b2c2e_b2mu_diffb/"); // LS
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/", "Pair/bbbar/b2c2mu_b2e_sameb/");
    fRegistry.addClone("Pair/ccbar/c2e_c2mu/", "Pair/bbbar/b2c2mu_b2e_diffb/"); // LS

    // for electron info
    fRegistry.add("Track/Electron/c2e/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/Electron/c2e/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/Electron/c2e/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/Electron/c2e/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/Electron/c2e/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/Electron/c2e/hDCAxy_Pt", "DCA_{xy} vs. pT;p_{T} (GeV/c);DCA_{xy} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/Electron/c2e/hDCAz_Pt", "DCA_{z} vs. pT;p_{T} (GeV/c);DCA_{z} (cm)", kTH2F, {{1000, 0, 10}, {200, -1, 1}}, false);
    fRegistry.add("Track/Electron/c2e/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/Electron/c2e/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/Electron/c2e/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/Electron/c2e/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/Electron/c2e/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTOFbeta", "TOF beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    fRegistry.add("Track/Electron/c2e/h1overTOFbeta", "TOF beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {1000, 0.8, 1.8}}, false);
    fRegistry.add("Track/Electron/c2e/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/Electron/c2e/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/Electron/c2e/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/Electron/c2e/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/Electron/c2e/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/Electron/c2e/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH1F, {{32, 0, 16}}, false);
    fRegistry.add("Track/Electron/c2e/hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/Electron/c2e/hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/Electron/c2e/hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.addClone("Track/Electron/c2e/", "Track/Electron/b2e/");
    fRegistry.addClone("Track/Electron/c2e/", "Track/Electron/b2c2e/");

    // for muon info
    fRegistry.add("Track/Muon/c2mu/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/Muon/c2mu/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/Muon/c2mu/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {25, -4.5f, -2.0f}}, false);
    fRegistry.add("Track/Muon/c2mu/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
    fRegistry.add("Track/Muon/c2mu/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/Muon/c2mu/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/Muon/c2mu/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/Muon/c2mu/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/Muon/c2mu/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
    fRegistry.add("Track/Muon/c2mu/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
    fRegistry.add("Track/Muon/c2mu/hPDCA", "pDCA;p_{T} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/Muon/c2mu/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/Muon/c2mu/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/Muon/c2mu/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/Muon/c2mu/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
    fRegistry.add("Track/Muon/c2mu/hPtGen_DeltaPtOverPtGen", "muon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/Muon/c2mu/hPtGen_DeltaEta", "muon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/Muon/c2mu/hPtGen_DeltaPhi", "muon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.addClone("Track/Muon/c2mu/", "Track/Muon/b2mu/");
    fRegistry.addClone("Track/Muon/c2mu/", "Track/Muon/b2c2mu/");

    // for generated info
    fRegistry.add("Generated/ccbar/c2e_c2mu/hadron_hadron/hs", "gen. e#mu correlation", kTHnSparseD, {axis_mass, axis_pt, axis_rapidity, axis_dphi}, true);
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/hadron_hadron/", "Generated/ccbar/c2e_c2mu/meson_meson/");
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/hadron_hadron/", "Generated/ccbar/c2e_c2mu/baryon_baryon/");
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/hadron_hadron/", "Generated/ccbar/c2e_c2mu/meson_baryon/");
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/", "Generated/bbbar/b2e_b2mu/");
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/", "Generated/bbbar/b2c2e_b2c2mu/");
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/", "Generated/bbbar/b2c2e_b2mu_sameb/");
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/", "Generated/bbbar/b2c2e_b2mu_diffb/"); // LS
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/", "Generated/bbbar/b2c2mu_b2e_sameb/");
    fRegistry.addClone("Generated/ccbar/c2e_c2mu/", "Generated/bbbar/b2c2mu_b2e_diffb/"); // LS
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

  template <typename TCollision, typename TElectron, typename TMuon, typename TMCParticles>
  bool fillTruePairInfo(TCollision const& collision, TElectron const& t1, TMuon const& t2, TMCParticles const& mcparticles)
  {
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

    auto t1mc = t1.template emmcparticle_as<TMCParticles>();
    auto t2mc = t2.template emmcparticle_as<TMCParticles>();
    if (abs(t1mc.pdgCode()) != 11 || abs(t2mc.pdgCode()) != 13) {
      return false;
    }

    if (!t1mc.isPhysicalPrimary() && !t1mc.producedByGenerator()) {
      return false;
    }
    if (!t2mc.isPhysicalPrimary() && !t2mc.producedByGenerator()) {
      return false;
    }

    if (t1mc.emmceventId() != t2mc.emmceventId()) {
      return false;
    }

    int hfee_type = IsHF(t1mc, t2mc, mcparticles);
    if (hfee_type < 0) {
      return false;
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

    auto mp1 = mcparticles.iteratorAt(t1mc.mothersIds()[0]);
    auto mp2 = mcparticles.iteratorAt(t2mc.mothersIds()[0]);
    if (t1mc.pdgCode() * t2mc.pdgCode() < 0) { // ULS
      switch (hfee_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce): {
          fRegistry.fill(HIST("Pair/ccbar/c2e_c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
          if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
            fRegistry.fill(HIST("Pair/ccbar/c2e_c2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<0, TMCParticles>(t1);
            fillMuonInfo<3, TMCParticles>(t2);
          } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
            fRegistry.fill(HIST("Pair/ccbar/c2e_c2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<0, TMCParticles>(t1);
            fillMuonInfo<3, TMCParticles>(t2);
          } else {
            fRegistry.fill(HIST("Pair/ccbar/c2e_c2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<0, TMCParticles>(t1);
            fillMuonInfo<3, TMCParticles>(t2);
          }
          break;
        }
        case static_cast<int>(EM_HFeeType::kBe_Be): {
          fRegistry.fill(HIST("Pair/bbbar/b2e_b2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
          if (isBeautyMeson(mp1) && isBeautyMeson(mp2)) {
            fRegistry.fill(HIST("Pair/bbbar/b2e_b2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<1, TMCParticles>(t1);
            fillMuonInfo<4, TMCParticles>(t2);
          } else if (isBeautyBaryon(mp1) && isBeautyBaryon(mp2)) {
            fRegistry.fill(HIST("Pair/bbbar/b2e_b2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<1, TMCParticles>(t1);
            fillMuonInfo<4, TMCParticles>(t2);
          } else {
            fRegistry.fill(HIST("Pair/bbbar/b2e_b2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<1, TMCParticles>(t1);
            fillMuonInfo<4, TMCParticles>(t2);
          }
          break;
        }
        case static_cast<int>(EM_HFeeType::kBCe_BCe): {
          fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
          if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
            fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<2, TMCParticles>(t1);
            fillMuonInfo<5, TMCParticles>(t2);
          } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
            fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<2, TMCParticles>(t1);
            fillMuonInfo<5, TMCParticles>(t2);
          } else {
            fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2c2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            fillElectronInfo<2, TMCParticles>(t1);
            fillMuonInfo<5, TMCParticles>(t2);
          }
          break;
        }
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): { // ULS
          if (isBeautyMeson(mp2) || isBeautyBaryon(mp2)) {
            fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_sameb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else {
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            }
            fillElectronInfo<2, TMCParticles>(t1);
            fillMuonInfo<4, TMCParticles>(t2);
          } else {
            fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_sameb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else {
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            }
            fillElectronInfo<1, TMCParticles>(t1);
            fillMuonInfo<5, TMCParticles>(t2);
          }
          break;
        }
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
          LOGF(info, "You should not see kBCe_Be_DiffB in ULS. Good luck.");
          break;
        default:
          break;
      }
    } else { // LS
      switch (hfee_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce):
          LOGF(info, "You should not see kCe_Ce in LS. Good luck.");
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be):
          LOGF(info, "You should not see kBe_Be in LS. Good luck.");
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe):
          LOGF(info, "You should not see kBCe_BCe in LS. Good luck.");
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
          LOGF(info, "You should not see kBCe_Be_SameB in LS. Good luck.");
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): { // LS
          if (isBeautyMeson(mp2) || isBeautyBaryon(mp2)) {
            fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_diffb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else {
              fRegistry.fill(HIST("Pair/bbbar/b2c2e_b2mu_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            }
            fillElectronInfo<2, TMCParticles>(t1);
            fillMuonInfo<4, TMCParticles>(t2);
          } else {
            fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_diffb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            } else {
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2e_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi, dca_3d_el, dca_xy_mu);
            }
            fillElectronInfo<1, TMCParticles>(t1);
            fillMuonInfo<5, TMCParticles>(t2);
          }
          break;
        }
        default:
          break;
      }
    }

    fillElectronInfo<0, TMCParticles>(t1);
    fillMuonInfo<0, TMCParticles>(t2);
    return true;
  }

  template <int e_source_id, typename TMCParticles, typename TTrack>
  void fillElectronInfo(TTrack const& track)
  {
    if (std::find(used_electronIds.begin(), used_electronIds.end(), track.globalIndex()) == used_electronIds.end()) {
      auto mctrack = track.template emmcparticle_as<TMCParticles>();
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hPt"), track.pt());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAxy_Pt"), track.pt(), track.dcaXY());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAz_Pt"), track.pt(), track.dcaZ());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hMeanClusterSizeITS"), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTOFbeta"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("h1overTOFbeta"), track.tpcInnerParam(), 1. / track.beta());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
      fRegistry.fill(HIST("Track/Electron/") + HIST(lepton_source_types[e_source_id]) + HIST("hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());
      used_electronIds.emplace_back(track.globalIndex());
    }
  }

  template <int e_source_id, typename TMCParticles, typename TTrack>
  void fillMuonInfo(TTrack const& track)
  {
    if (std::find(used_muonIds.begin(), used_muonIds.end(), track.globalIndex()) == used_muonIds.end()) {
      auto mctrack = track.template emmcparticle_as<TMCParticles>();
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hPt"), track.pt());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hTrackType"), track.trackType());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAxySigma"), track.fwdDcaX() / sqrt(track.cXX()), track.fwdDcaY() / sqrt(track.cYY()));
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAxRes_Pt"), track.pt(), sqrt(track.cXX()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hDCAyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hNclsMFT"), track.nClustersMFT());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hNclsMCH"), track.nClusters());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hPDCA"), track.pt(), track.pDca());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hChi2"), track.chi2());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hChi2MatchMCHMID"), track.chi2MatchMCHMID());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hMFTClusterMap"), track.mftClusterMap());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
      fRegistry.fill(HIST("Track/Muon/") + HIST(lepton_source_types[e_source_id]) + HIST("hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());
      used_muonIds.emplace_back(track.globalIndex());
    }
  }

  template <typename T>
  bool isElectronInAcceptance(T const& t1)
  {
    if ((mctrackcuts.min_mcPt_electron < t1.pt() && t1.pt() < mctrackcuts.max_mcPt_electron) && mctrackcuts.min_mcEta_electron < t1.eta() && t1.eta() < mctrackcuts.max_mcEta_electron) {
      return true;
    } else {
      return false;
    }
  }

  template <typename T>
  bool isMuonInAcceptance(T const& t1)
  {
    if ((mctrackcuts.min_mcPt_muon < t1.pt() && t1.pt() < mctrackcuts.max_mcPt_muon) && mctrackcuts.min_mcEta_muon < t1.eta() && t1.eta() < mctrackcuts.max_mcEta_muon) {
      return true;
    } else {
      return false;
    }
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

  Partition<FilteredMyElectrons> posTracks_el = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyElectrons> negTracks_el = o2::aod::emprimaryelectron::sign < int8_t(0);
  Partition<FilteredMyMuons> posTracks_mu = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<FilteredMyMuons> negTracks_mu = o2::aod::emprimarymuon::sign < int8_t(0);

  std::vector<int> used_electronIds;
  std::vector<int> used_muonIds;
  void processQCMC(FilteredMyCollisions const& collisions, FilteredMyElectrons const&, FilteredMyMuons const&, aod::EMMCParticles const& mcparticles)
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

      for (auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_el_per_coll, negTracks_mu_per_coll))) { // ULS+-
        fillTruePairInfo(collision, pos, neg, mcparticles);
      }

      for (auto& [neg, pos] : combinations(CombinationsFullIndexPolicy(negTracks_el_per_coll, posTracks_mu_per_coll))) { // ULS-+
        fillTruePairInfo(collision, neg, pos, mcparticles);
      }

      for (auto& [pos1, pos2] : combinations(CombinationsFullIndexPolicy(posTracks_el_per_coll, posTracks_mu_per_coll))) { // LS++
        fillTruePairInfo(collision, pos1, pos2, mcparticles);
      }
      for (auto& [neg1, neg2] : combinations(CombinationsFullIndexPolicy(negTracks_el_per_coll, negTracks_mu_per_coll))) { // LS--
        fillTruePairInfo(collision, neg1, neg2, mcparticles);
      }

    } // end of collision loop

    used_electronIds.clear();
    used_electronIds.shrink_to_fit();
    used_muonIds.clear();
    used_muonIds.shrink_to_fit();
  } // end of process
  PROCESS_SWITCH(emuCorrelationMC, processQCMC, "run dielectron QCMC", true);

  Partition<aod::EMMCParticles> pos_els_MC = o2::aod::mcparticle::pdgCode == -11; // e+
  Partition<aod::EMMCParticles> neg_els_MC = o2::aod::mcparticle::pdgCode == +11; // e-
  Partition<aod::EMMCParticles> pos_mus_MC = o2::aod::mcparticle::pdgCode == -13; // mu+
  Partition<aod::EMMCParticles> neg_mus_MC = o2::aod::mcparticle::pdgCode == +13; // mu-
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(MyCollisions const& collisions, aod::EMMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      auto mccollision = collision.emmcevent_as<aod::EMMCEvents>();
      auto pos_els_per_coll = pos_els_MC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
      auto neg_els_per_coll = neg_els_MC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
      auto pos_mus_per_coll = pos_mus_MC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
      auto neg_mus_per_coll = neg_mus_MC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);

      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(pos_els_per_coll, neg_mus_per_coll))) { // ULS+-

        if (!isElectronInAcceptance(t1) || !isMuonInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int hfee_type = IsHF(t1, t2, mcparticles);
        if (hfee_type < 0) {
          continue;
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

        auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
        auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
        switch (hfee_type) {
          case static_cast<int>(EM_HFeeType::kCe_Ce): {
            fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
              fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
              fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else {
              fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBe_Be): {
            fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            if (isBeautyMeson(mp1) && isBeautyMeson(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else if (isBeautyBaryon(mp1) && isBeautyBaryon(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBCe_BCe): {
            fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): { // ULS
            if (isBeautyMeson(mp2) || isBeautyBaryon(mp2)) {   // muon is from b hadron.
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            } else { // electron is from b hadron
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
            LOGF(info, "You should not see kBCe_Be_DiffB in ULS. Good luck.");
            break;
          default:
            break;
        }
      } // end of true ULS pair loop

      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(neg_els_per_coll, pos_mus_per_coll))) { // ULS-+

        if (!isElectronInAcceptance(t1) || !isMuonInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int hfee_type = IsHF(t1, t2, mcparticles);
        if (hfee_type < 0) {
          continue;
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

        auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
        auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
        switch (hfee_type) {
          case static_cast<int>(EM_HFeeType::kCe_Ce): {
            fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
              fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
              fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else {
              fRegistry.fill(HIST("Generated/ccbar/c2e_c2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBe_Be): {
            fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            if (isBeautyMeson(mp1) && isBeautyMeson(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else if (isBeautyBaryon(mp1) && isBeautyBaryon(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBCe_BCe): {
            fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            } else {
              fRegistry.fill(HIST("Generated/bbbar/b2e_b2mu/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): { // ULS
            if (isBeautyMeson(mp2) || isBeautyBaryon(mp2)) {   // muon is from b hadron.
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            } else { // electron is from b hadron
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            }
            break;
          }
          case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
            LOGF(info, "You should not see kBCe_Be_DiffB in ULS. Good luck.");
            break;
          default:
            break;
        }
      } // end of true ULS pair loop

      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(pos_els_per_coll, pos_mus_per_coll))) {
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isElectronInAcceptance(t1) || !isMuonInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int hfee_type = IsHF(t1, t2, mcparticles);
        if (hfee_type < 0) {
          continue;
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
        auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
        auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
        switch (hfee_type) {
          case static_cast<int>(EM_HFeeType::kCe_Ce):
            LOGF(info, "You should not see kCe_Ce in LS++. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBe_Be):
            LOGF(info, "You should not see kBe_Be in LS++. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBCe_BCe):
            LOGF(info, "You should not see kBCe_BCe in LS++. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
            LOGF(info, "You should not see kBCe_Be_SameB in LS++. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): { // LS
            if (isBeautyMeson(mp2) || isBeautyBaryon(mp2)) {   // muon is from b hadron.
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            } else { // electron is from b hadron
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            }
            break;
          }
          default:
            break;
        }
      } // end of true LS++ pair loop

      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(neg_els_per_coll, neg_mus_per_coll))) {
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isElectronInAcceptance(t1) || !isMuonInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int hfee_type = IsHF(t1, t2, mcparticles);
        if (hfee_type < 0) {
          continue;
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
        auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
        auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
        switch (hfee_type) {
          case static_cast<int>(EM_HFeeType::kCe_Ce):
            LOGF(info, "You should not see kCe_Ce in LS--. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBe_Be):
            LOGF(info, "You should not see kBe_Be in LS--. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBCe_BCe):
            LOGF(info, "You should not see kBCe_BCe in LS--. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
            LOGF(info, "You should not see kBCe_Be_SameB in LS--. Good luck.");
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): { // LS
            if (isBeautyMeson(mp2) || isBeautyBaryon(mp2)) {   // muon is from b hadron.
              fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2e_b2mu_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            } else { // electron is from b hadron
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/meson_meson/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2e_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), v12.Rapidity(), dphi);
              }
            }
            break;
          }
          default:
            break;
        }
      } // end of true LS++ pair loop
    }   // end of collision loop
  }
  PROCESS_SWITCH(emuCorrelationMC, processGen, "run genrated info", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(emuCorrelationMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<emuCorrelationMC>(cfgc, TaskName{"emu-correlation-mc"})};
}
