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
// This code runs loop over electrons for QC in MC.
//    Please write to: daiki.sekihata@cern.ch

#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit, aod::EMPrimaryElectronMCLabels>;
using MyMCTrack = MyMCTracks::iterator;

struct singleElectronQCMC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  ConfigurableAxis ConfPteBins{"ConfPteBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00}, "pTe bins for output histograms"};

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

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.8, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};

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
    Configurable<bool> enableTTCA{"enableTTCA", false, "Flag to enable or disable TTCA"}; // TTCA should not be used for single track to avoid double counting.

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

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt{"min_mcPt", 0.05, "min. MC pT"};
    Configurable<float> max_mcPt{"max_mcPt", 1e+10, "max. MC pT"};
    Configurable<float> min_mcEta{"min_mcEta", -0.9, "max. MC eta"};
    Configurable<float> max_mcEta{"max_mcEta", +0.9, "max. MC eta"};
  } mctrackcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false}; // 1 HistogramRegistry can keep up to 512 histograms
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view ele_source_types[9] = {"lf/", "Photon/", "PromptJPsi/", "NonPromptJPsi/", "PromptPsi2S/", "NonPromptPsi2S/", "c2e/", "b2e/", "b2c2e/"};

  ~singleElectronQCMC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);

    const AxisSpec axis_pt{ConfPteBins, "p_{T,e} (GeV/c)"};
    const AxisSpec axis_eta{20, -1.0, +1.0, "#eta_{e}"};
    const AxisSpec axis_phi{36, 0.0, 2 * M_PI, "#varphi_{e} (rad.)"};
    const AxisSpec axis_charge_rec{3, -1.5, +1.5, "charge"};
    const AxisSpec axis_charge_gen{3, -1.5, +1.5, "true charge"};

    // generated info
    fRegistry.add("Generated/lf/hs", "gen. single electron", kTHnSparseF, {axis_pt, axis_eta, axis_phi, axis_charge_gen}, true);
    fRegistry.addClone("Generated/lf/", "Generated/PromptJPsi/");
    fRegistry.addClone("Generated/lf/", "Generated/NonPromptJPsi/");
    fRegistry.addClone("Generated/lf/", "Generated/PromptPsi2S/");
    fRegistry.addClone("Generated/lf/", "Generated/NonPromptPsi2S/");
    fRegistry.addClone("Generated/lf/", "Generated/c2e/");
    fRegistry.addClone("Generated/lf/", "Generated/b2e/");
    fRegistry.addClone("Generated/lf/", "Generated/b2c2e/");

    // track info
    fRegistry.add("Track/lf/hs", "rec. single electron", kTHnSparseF, {axis_pt, axis_eta, axis_phi, axis_charge_rec, axis_charge_gen}, true);
    fRegistry.add("Track/lf/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/lf/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/lf/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/lf/hDCA3DSigma", "DCA 3D;DCA_{3D} (#sigma);", kTH1F, {{100, 0.0f, 10.0f}}, false);
    fRegistry.add("Track/lf/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/lf/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/lf/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/lf/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/lf/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/lf/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFbeta", "TOF #beta;p_{in} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    fRegistry.add("Track/lf/h1overTOFbeta", "TOF 1/#beta;p_{in} (GeV/c);1/#beta", kTH2F, {{1000, 0, 10}, {1000, 0.8, 1.8}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/lf/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/lf/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/lf/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/lf/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/lf/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/lf/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH1F, {{32, 0, 16}}, false);
    fRegistry.add("Track/lf/hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/lf/hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/lf/hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.addClone("Track/lf/", "Track/Photon/"); // this is not for efficiency! only for contamination. We don't store generated photon conversions.
    fRegistry.addClone("Track/lf/", "Track/PromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/PromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/c2e/");
    fRegistry.addClone("Track/lf/", "Track/b2e/");
    fRegistry.addClone("Track/lf/", "Track/b2c2e/");
  }

  bool cfgDoFlow = false;
  void init(InitContext&)
  {
    DefineEMEventCut();
    DefineDileptonCut();
    addhistograms();

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
    fDielectonCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectonCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);

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

  template <typename T>
  bool isInAcceptance(T const& t1)
  {
    if ((mctrackcuts.min_mcPt < t1.pt() && t1.pt() < mctrackcuts.max_mcPt) && (mctrackcuts.min_mcEta < t1.eta() && t1.eta() < mctrackcuts.max_mcEta)) {
      return true;
    } else {
      return false;
    }
  }

  template <int e_source_id, typename TMCParticles, typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    auto mctrack = track.template emmcparticle_as<TMCParticles>();
    float dca_3d = dca3DinSigma(track);
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hs"), track.pt(), track.eta(), track.phi(), track.sign(), -mctrack.pdgCode() / 11);
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxyz"), track.dcaXY(), track.dcaZ());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCA3DSigma"), dca_3d);
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNclsITS"), track.itsNCls());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNclsTPC"), track.tpcNClsFound());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNcrTPC"), track.tpcNClsCrossedRows());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hChi2TPC"), track.tpcChi2NCl());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hChi2ITS"), track.itsChi2NCl());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hITSClusterMap"), track.itsClusterMap());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hMeanClusterSizeITS"), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFbeta"), track.tpcInnerParam(), track.beta());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("h1overTOFbeta"), track.tpcInnerParam(), 1. / track.beta());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
    fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());
  }

  std::vector<int> used_trackIds;
  SliceCache cache;
  Preslice<MyMCTracks> perCollision_track = aod::emprimaryelectron::emeventId;
  Filter trackFilter = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dielectroncuts.cfg_max_eta_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter = (dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < dielectroncuts.cfg_min_TPCNsigmaPi || dielectroncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi) && ((0.96f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.04f) || o2::aod::pidtofbeta::beta < 0.f);
  Filter ttcaFilter = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  using FilteredMyMCTracks = soa::Filtered<MyMCTracks>;

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processQCMC(FilteredMyCollisions const& collisions, FilteredMyMCTracks const& tracks, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const&)
  {
    used_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      initCCDB(collision);
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision, cfgDoFlow);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision, cfgDoFlow);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto mccollision = collision.emmcevent_as<aod::EMMCEvents>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCollision_track, collision.globalIndex());

      for (auto& track : tracks_per_coll) {
        auto mctrack = track.template emmcparticle_as<aod::EMMCParticles>();
        if (abs(mctrack.pdgCode()) != 11) {
          continue;
        }
        if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
          if (!fDielectonCut.IsSelectedTrack<true>(track, collision)) {
            continue;
          }
        } else { // cut-based
          if (!fDielectonCut.IsSelectedTrack(track)) {
            continue;
          }
        }

        if (!mctrack.has_mothers()) {
          continue;
        }
        auto mcmother = mcparticles.iteratorAt(mctrack.mothersIds()[0]);
        int pdg_mother = abs(mcmother.pdgCode());

        if (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) {
          if (pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 331 || pdg_mother == 113 || pdg_mother == 223 || pdg_mother == 333) {
            fillTrackInfo<0, aod::EMMCParticles>(track);
          } else if (pdg_mother == 443) {
            if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
              fillTrackInfo<3, aod::EMMCParticles>(track);
            } else {
              fillTrackInfo<2, aod::EMMCParticles>(track);
            }
          } else if (pdg_mother == 100443) {
            if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
              fillTrackInfo<5, aod::EMMCParticles>(track);
            } else {
              fillTrackInfo<4, aod::EMMCParticles>(track);
            }
          } else if (IsFromBeauty(mctrack, mcparticles) > 0) { // b is found in full decay chain.
            if (IsFromCharm(mctrack, mcparticles) > 0) {       // c is found in full decay chain.
              fillTrackInfo<8, aod::EMMCParticles>(track);
            } else {
              fillTrackInfo<7, aod::EMMCParticles>(track);
            }
          } else if (IsFromCharm(mctrack, mcparticles) > 0) { // c is found in full decay chain. Not from b.
            fillTrackInfo<6, aod::EMMCParticles>(track);
          }
        } else {
          fillTrackInfo<1, aod::EMMCParticles>(track);
        }

      } // end of track loop

    } // end of collision loop

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();

  } // end of process
  PROCESS_SWITCH(singleElectronQCMC, processQCMC, "run dielectron QC MC", true);

  Partition<aod::EMMCParticles> electronsMC = nabs(o2::aod::mcparticle::pdgCode) == 11; // e+, e-
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
      // LOGF(info, "mccollision.getGeneratorId() = %d", mccollision.getGeneratorId());
      // LOGF(info, "mccollision.getSubGeneratorId() = %d", mccollision.getSubGeneratorId());
      // LOGF(info, "mccollision.getSourceId() = %d", mccollision.getSourceId());
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      auto electronsMC_per_coll = electronsMC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);

      for (auto& ele : electronsMC_per_coll) {
        if (!(ele.isPhysicalPrimary() || ele.producedByGenerator())) {
          continue;
        }
        if (!isInAcceptance(ele)) {
          continue;
        }
        if (!ele.has_mothers()) {
          continue;
        }
        auto mcmother = mcparticles.iteratorAt(ele.mothersIds()[0]);
        int pdg_mother = abs(mcmother.pdgCode());
        if (pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 331 || pdg_mother == 113 || pdg_mother == 223 || pdg_mother == 333) {
          fRegistry.fill(HIST("Generated/lf/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
        } else if (pdg_mother == 443) {
          if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
            fRegistry.fill(HIST("Generated/NonPromptJPsi/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
          } else {
            fRegistry.fill(HIST("Generated/PromptJPsi/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
          }
        } else if (pdg_mother == 100443) {
          if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
            fRegistry.fill(HIST("Generated/NonPromptPsi2S/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
          } else {
            fRegistry.fill(HIST("Generated/PromptPsi2S/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
          }
        } else if (IsFromBeauty(ele, mcparticles) > 0) { // b is found in full decay chain.
          if (IsFromCharm(ele, mcparticles) > 0) {       // c is found in full decay chain.
            fRegistry.fill(HIST("Generated/b2c2e/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
          } else {
            fRegistry.fill(HIST("Generated/b2e/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
          }
        } else if (IsFromCharm(ele, mcparticles) > 0) { // c is found in full decay chain. Not from b.
          fRegistry.fill(HIST("Generated/c2e/hs"), ele.pt(), ele.eta(), ele.phi(), -ele.pdgCode() / 11);
        }
      }

    } // end of collision loop
  }
  PROCESS_SWITCH(singleElectronQCMC, processGen, "run genrated info", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(singleElectronQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<singleElectronQCMC>(cfgc, TaskName{"single-electron-qc-mc"})};
}
