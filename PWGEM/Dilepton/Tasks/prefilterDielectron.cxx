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
// This code produces information on prefilter for dielectron.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TString.h"

#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;

struct prefilterDielectron {
  Produces<aod::EMPrimaryElectronsPrefilterBitDerived> pfb_derived;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
  } eventcuts;

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "dielectroncut_group";

    // for mee prefilter
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass for prefilter ULS"}; // region to be rejected
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.0, "max mass for prefilter ULS"}; // region to be rejected

    // for phiv prefilter
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", false, "flag to apply phiv cut"};              // region to be rejected
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};              // region to be rejected
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"}; // region to be rejected
    Configurable<float> cfg_min_phiv{"cfg_min_phiv", -1.f, "min phiv"};                                // region to be rejected
    Configurable<float> cfg_max_phiv{"cfg_max_phiv", 3.2, "max phiv"};                                 // region to be rejected

    // for deta-dphi prefilter
    Configurable<bool> cfg_apply_detadphi_uls{"cfg_apply_detadphi_uls", false, "flag to apply generator deta-dphi elliptic cut in ULS"}; // region to be rejected
    Configurable<bool> cfg_apply_detadphi_ls{"cfg_apply_detadphi_ls", false, "flag to apply generator deta-dphi elliptic cut in LS"};    // region to be rejected
    Configurable<float> cfg_min_deta_ls{"cfg_min_deta_ls", 0.04, "deta between 2 electrons (elliptic cut)"};                             // region to be rejected
    Configurable<float> cfg_min_dphi_ls{"cfg_min_dphi_ls", 0.2, "dphi between 2 electrons (elliptic cut)"};                              // region to be rejected
    Configurable<float> cfg_min_deta_uls{"cfg_min_deta_uls", 0.04, "deta between 2 electrons (elliptic cut)"};                           // region to be rejected
    Configurable<float> cfg_min_dphi_uls{"cfg_min_dphi_uls", 0.2, "dphi between 2 electrons (elliptic cut)"};                            // region to be rejected

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.15, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "min phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2 TOF"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.f, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.f, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_min_its_cluster_size{"cfg_min_its_cluster_size", 0.f, "min ITS cluster size"};
    Configurable<float> cfg_max_its_cluster_size{"cfg_max_its_cluster_size", 16.f, "max ITS cluster size"};
    Configurable<float> cfg_min_rel_diff_pin{"cfg_min_rel_diff_pin", -1e+10, "min rel. diff. between pin and ppv"};
    Configurable<float> cfg_max_rel_diff_pin{"cfg_max_rel_diff_pin", +1e+10, "max rel. diff. between pin and ppv"};
    // Configurable<float> cfgRefR{"cfgRefR", 1.2, "reference R (in m) for extrapolation"}; // https://cds.cern.ch/record/1419204

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3, kTOFif : 4, kPIDML : 5, kTPChadrejORTOFreq_woTOFif : 6]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    // Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    // Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
    Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to enable ITSsa tracks"};
    Configurable<float> cfg_max_pt_track_ITSsa{"cfg_max_pt_track_ITSsa", 0.15, "max pt for ITSsa tracks"};

    // configuration for PID ML
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
    Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
    Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{-999999., 999999.}, "Bin limits for ML application"};
    Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.95}, "ML cuts per bin"};
    Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature"}, "Names of ML model input features"};
    Configurable<std::string> nameBinningFeature{"nameBinningFeature", "pt", "Names of ML model binning feature"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dielectroncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext& /*context*/)
  {
    DefineEMEventCut();
    DefineDielectronCut();
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
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = collision.runNumber();
  }

  ~prefilterDielectron() {}

  void addhistograms()
  {
    const AxisSpec axis_mass{400, 0, 4, "m_{ee} (GeV/c^{2})"};
    const AxisSpec axis_pair_pt{100, 0, 10, "p_{T,ee} (GeV/c)"};
    const AxisSpec axis_phiv{90, 0, M_PI, "#varphi_{V} (rad.)"};

    // for pair
    fRegistry.add("Pair/before/uls/hMvsPt", "m_{ee} vs. p_{T,ee}", kTH2D, {axis_mass, axis_pair_pt}, true);
    fRegistry.add("Pair/before/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2D, {axis_phiv, {200, 0, 1}}, true);
    fRegistry.add("Pair/before/uls/hDeltaEtaDeltaPhi", "#Delta#eta-#Delta#varphi between 2 tracks;#Delta#varphi (rad.);#Delta#eta;", kTH2D, {{180, -M_PI, M_PI}, {400, -2, +2}}, true);
    fRegistry.addClone("Pair/before/uls/", "Pair/before/lspp/");
    fRegistry.addClone("Pair/before/uls/", "Pair/before/lsmm/");
    fRegistry.addClone("Pair/before/", "Pair/after/");
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
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
  }

  o2::analysis::MlResponseDielectronSingleTrack<float> mlResponseSingleTrack;
  void DefineDielectronCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for pair
    fDielectronCut.SetMeeRange(0, 1e+10);
    fDielectronCut.SetPairPtRange(0.f, 1e+10);
    fDielectronCut.SetPairYRange(-1e+10, +1e+10);
    fDielectronCut.SetPairDCARange(0.f, 1e+10); // in sigma
    fDielectronCut.ApplyPhiV(false);
    fDielectronCut.ApplyPrefilter(false);
    fDielectronCut.SetMindEtadPhi(false, false, 1.f, 1.f);
    fDielectronCut.SetPairOpAng(0.f, 3.2f);
    fDielectronCut.SetRequireDifferentSides(false);

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, dielectroncuts.cfg_max_pt_track);
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetTrackPhiRange(dielectroncuts.cfg_min_phi_track, dielectroncuts.cfg_max_phi_track);
    fDielectronCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetMaxFracSharedClustersTPC(dielectroncuts.cfg_max_frac_shared_clusters_tpc);
    fDielectronCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITS(dielectroncuts.cfg_min_its_cluster_size, dielectroncuts.cfg_max_its_cluster_size);
    fDielectronCut.SetTrackMaxDcaXY(dielectroncuts.cfg_max_dcaxy);
    fDielectronCut.SetTrackMaxDcaZ(dielectroncuts.cfg_max_dcaz);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);
    fDielectronCut.SetChi2TOF(0, dielectroncuts.cfg_max_chi2tof);
    fDielectronCut.SetRelDiffPin(dielectroncuts.cfg_min_rel_diff_pin, dielectroncuts.cfg_max_rel_diff_pin);
    fDielectronCut.IncludeITSsa(dielectroncuts.includeITSsa, dielectroncuts.cfg_max_pt_track_ITSsa);

    // for eID
    fDielectronCut.SetPIDScheme(dielectroncuts.cfg_pid_scheme);
    fDielectronCut.SetTPCNsigmaElRange(dielectroncuts.cfg_min_TPCNsigmaEl, dielectroncuts.cfg_max_TPCNsigmaEl);
    // fDielectronCut.SetTPCNsigmaMuRange(dielectroncuts.cfg_min_TPCNsigmaMu, dielectroncuts.cfg_max_TPCNsigmaMu);
    fDielectronCut.SetTPCNsigmaPiRange(dielectroncuts.cfg_min_TPCNsigmaPi, dielectroncuts.cfg_max_TPCNsigmaPi);
    fDielectronCut.SetTPCNsigmaKaRange(dielectroncuts.cfg_min_TPCNsigmaKa, dielectroncuts.cfg_max_TPCNsigmaKa);
    fDielectronCut.SetTPCNsigmaPrRange(dielectroncuts.cfg_min_TPCNsigmaPr, dielectroncuts.cfg_max_TPCNsigmaPr);
    fDielectronCut.SetTOFNsigmaElRange(dielectroncuts.cfg_min_TOFNsigmaEl, dielectroncuts.cfg_max_TOFNsigmaEl);

    if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      std::vector<float> binsML{};
      binsML.reserve(dielectroncuts.binsMl.value.size());
      for (size_t i = 0; i < dielectroncuts.binsMl.value.size(); i++) {
        binsML.emplace_back(dielectroncuts.binsMl.value[i]);
      }
      std::vector<float> thresholdsML{};
      thresholdsML.reserve(dielectroncuts.cutsMl.value.size());
      for (size_t i = 0; i < dielectroncuts.cutsMl.value.size(); i++) {
        thresholdsML.emplace_back(dielectroncuts.cutsMl.value[i]);
      }
      fDielectronCut.SetMLThresholds(binsML, thresholdsML);

      // static constexpr int nClassesMl = 2;
      // const std::vector<int> cutDirMl = {o2::cuts_ml::CutSmaller, o2::cuts_ml::CutNot};
      // const std::vector<std::string> labelsClasses = {"Signal", "Background"};
      // const uint32_t nBinsMl = dielectroncuts.binsMl.value.size() - 1;
      // const std::vector<std::string> labelsBins(nBinsMl, "bin");
      // double cutsMlArr[nBinsMl][nClassesMl];
      // for (uint32_t i = 0; i < nBinsMl; i++) {
      //   cutsMlArr[i][0] = dielectroncuts.cutsMl.value[i];
      //   cutsMlArr[i][1] = 0.;
      // }
      // o2::framework::LabeledArray<double> cutsMl = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      // mlResponseSingleTrack.configure(dielectroncuts.binsMl.value, cutsMl, cutDirMl, nClassesMl);
      // if (dielectroncuts.loadModelsFromCCDB) {
      //   ccdbApi.init(ccdburl);
      //   mlResponseSingleTrack.setModelPathsCCDB(dielectroncuts.onnxFileNames.value, ccdbApi, dielectroncuts.onnxPathsCCDB.value, dielectroncuts.timestampCCDB.value);
      // } else {
      //   mlResponseSingleTrack.setModelPathsLocal(dielectroncuts.onnxFileNames.value);
      // }
      // mlResponseSingleTrack.cacheInputFeaturesIndices(dielectroncuts.namesInputFeatures);
      // mlResponseSingleTrack.cacheBinningIndex(dielectroncuts.nameBinningFeature);
      // mlResponseSingleTrack.init(dielectroncuts.enableOptimizations.value);

      // fDielectronCut.SetPIDMlResponse(&mlResponseSingleTrack);
    } // end of PID ML
  }

  std::unordered_map<int, uint16_t> map_pfb; // map track.globalIndex -> prefilter bit

  SliceCache cache;
  Preslice<MyTracks> perCollision_track = aod::emprimaryelectron::emeventId;
  Partition<MyTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<MyTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  int ndf = 0;
  void processPFB(FilteredMyCollisions const& collisions, MyTracks const& tracks)
  {
    for (auto& track : tracks) {
      map_pfb[track.globalIndex()] = 0;
    } // end of track loop

    for (auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      bool is_cent_ok = true;
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        is_cent_ok = false;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

      if (!fEMEventCut.IsSelected(collision) || !is_cent_ok) {
        for (auto& pos : posTracks_per_coll) {
          map_pfb[pos.globalIndex()] = 0;
        }
        for (auto& neg : negTracks_per_coll) {
          map_pfb[neg.globalIndex()] = 0;
        }
        continue;
      }

      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (!fDielectronCut.IsSelectedTrack(pos) || !fDielectronCut.IsSelectedTrack(ele)) {
          continue;
        }
        // don't apply pair cut when you produce prefilter bit.

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), pos.sign(), ele.sign(), d_bz);
        float deta = pos.sign() * v1.Pt() > ele.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos.sign() * v1.Pt() > ele.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/before/uls/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/before/uls/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/before/uls/hDeltaEtaDeltaPhi"), dphi, deta);

        if (dielectroncuts.cfg_min_mass < v12.M() && v12.M() < dielectroncuts.cfg_max_mass) {
          map_pfb[pos.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kMee);
          map_pfb[ele.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kMee);
        }

        if (dielectroncuts.cfg_apply_phiv && ((v12.M() < dielectroncuts.cfg_phiv_slope * phiv + dielectroncuts.cfg_phiv_intercept) && (dielectroncuts.cfg_min_phiv < phiv && phiv < dielectroncuts.cfg_max_phiv))) {
          map_pfb[pos.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kPhiV);
          map_pfb[ele.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kPhiV);
        }

        if (dielectroncuts.cfg_apply_detadphi_uls && std::pow(deta / dielectroncuts.cfg_min_deta_uls, 2) + std::pow(dphi / dielectroncuts.cfg_min_dphi_uls, 2) < 1.f) {
          map_pfb[pos.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS);
          map_pfb[ele.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackULS);
        }
      } // end of ULS pairing

      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (!fDielectronCut.IsSelectedTrack(pos1) || !fDielectronCut.IsSelectedTrack(pos2)) {
          continue;
        }
        // don't apply pair cut when you produce prefilter bit.

        ROOT::Math::PtEtaPhiMVector v1(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), pos1.sign(), pos2.sign(), d_bz);
        float deta = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/before/lspp/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/before/lspp/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/before/lspp/hDeltaEtaDeltaPhi"), dphi, deta);

        if (dielectroncuts.cfg_apply_detadphi_ls && std::pow(deta / dielectroncuts.cfg_min_deta_ls, 2) + std::pow(dphi / dielectroncuts.cfg_min_dphi_ls, 2) < 1.f) {
          map_pfb[pos1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
          map_pfb[pos2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
        }
      } // end of LS++ pairing

      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (!fDielectronCut.IsSelectedTrack(ele1) || !fDielectronCut.IsSelectedTrack(ele2)) {
          continue;
        }
        // don't apply pair cut when you produce prefilter bit.

        ROOT::Math::PtEtaPhiMVector v1(ele1.pt(), ele1.eta(), ele1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), ele1.sign(), ele2.sign(), d_bz);
        float deta = ele1.sign() * v1.Pt() > ele2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = ele1.sign() * v1.Pt() > ele2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/before/lsmm/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/before/lsmm/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/before/lsmm/hDeltaEtaDeltaPhi"), dphi, deta);

        if (dielectroncuts.cfg_apply_detadphi_ls && std::pow(deta / dielectroncuts.cfg_min_deta_ls, 2) + std::pow(dphi / dielectroncuts.cfg_min_dphi_ls, 2) < 1.f) {
          map_pfb[ele1.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
          map_pfb[ele2.globalIndex()] |= 1 << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBitDerived::kSplitOrMergedTrackLS);
        }
      } // end of LS-- pairing

    } // end of collision loop

    for (auto& track : tracks) {
      // LOGF(info, "map_pfb[%d] = %d", track.globalIndex(), map_pfb[track.globalIndex()]);
      pfb_derived(map_pfb[track.globalIndex()]);
    } // end of track loop

    // check pfb.
    for (auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        if (!fDielectronCut.IsSelectedTrack(pos) || !fDielectronCut.IsSelectedTrack(ele)) {
          continue;
        }
        if (map_pfb[pos.globalIndex()] != 0 || map_pfb[ele.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), pos.sign(), ele.sign(), d_bz);
        float deta = pos.sign() * v1.Pt() > ele.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos.sign() * v1.Pt() > ele.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/after/uls/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/after/uls/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/after/uls/hDeltaEtaDeltaPhi"), dphi, deta);
      }

      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        if (!fDielectronCut.IsSelectedTrack(pos1) || !fDielectronCut.IsSelectedTrack(pos2)) {
          continue;
        }
        if (map_pfb[pos1.globalIndex()] != 0 || map_pfb[pos2.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), pos1.sign(), pos2.sign(), d_bz);
        float deta = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = pos1.sign() * v1.Pt() > pos2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/after/lspp/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/after/lspp/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/after/lspp/hDeltaEtaDeltaPhi"), dphi, deta);
      }

      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        if (!fDielectronCut.IsSelectedTrack(ele1) || !fDielectronCut.IsSelectedTrack(ele2)) {
          continue;
        }
        if (map_pfb[ele1.globalIndex()] != 0 || map_pfb[ele2.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(ele1.pt(), ele1.eta(), ele1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), ele1.sign(), ele2.sign(), d_bz);
        float deta = ele1.sign() * v1.Pt() > ele2.sign() * v2.Pt() ? v1.Eta() - v2.Eta() : v2.Eta() - v1.Eta();
        float dphi = ele1.sign() * v1.Pt() > ele2.sign() * v2.Pt() ? v1.Phi() - v2.Phi() : v2.Phi() - v1.Phi();
        o2::math_utils::bringToPMPi(dphi);

        fRegistry.fill(HIST("Pair/after/lsmm/hMvsPt"), v12.M(), v12.Pt());
        fRegistry.fill(HIST("Pair/after/lsmm/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/after/lsmm/hDeltaEtaDeltaPhi"), dphi, deta);
      }

    } // end of collision loop
    ndf++;
    map_pfb.clear();
  } // end of process
  PROCESS_SWITCH(prefilterDielectron, processPFB, "produce prefilter bit", false);

  void processDummy(MyTracks const& tracks)
  {
    for (int i = 0; i < tracks.size(); i++) {
      pfb_derived(0);
    }
  }
  PROCESS_SWITCH(prefilterDielectron, processDummy, "dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<prefilterDielectron>(cfgc, TaskName{"prefilter-dielectron"})};
}
