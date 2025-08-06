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
// This code runs loop over ULS ee pars for virtual photon QC.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Core/DielectronCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrack.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/MlResponseDielectronSingleTrack.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/Core/RecoDecay.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TString.h"

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;

struct vpPairQC {
  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  ConfigurableAxis ConfPtlBins{"ConfPtlBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTl bins for output histograms"};

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
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
  } eventcuts;

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "dielectroncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.01, "max mass"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pT"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pT"};
    Configurable<float> cfg_min_pair_y{"cfg_min_pair_y", -0.9, "min pair rapidity"};
    Configurable<float> cfg_max_pair_y{"cfg_max_pair_y", +0.9, "max pair rapidity"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};
    Configurable<float> cfg_min_phiv{"cfg_min_phiv", 0.0, "min phiv (constant)"};
    Configurable<float> cfg_max_phiv{"cfg_max_phiv", 3.2, "max phiv (constant)"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
    Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1e+10, "max DCA 3D in sigma"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2 TOF"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<float> cfg_min_its_cluster_size{"cfg_min_its_cluster_size", 0.f, "min ITS cluster size"};
    Configurable<float> cfg_max_its_cluster_size{"cfg_max_its_cluster_size", 16.f, "max ITS cluster size"};
    Configurable<float> cfg_min_rel_diff_pin{"cfg_min_rel_diff_pin", -1e+10, "min rel. diff. between pin and ppv"};
    Configurable<float> cfg_max_rel_diff_pin{"cfg_max_rel_diff_pin", +1e+10, "max rel. diff. between pin and ppv"};
    Configurable<float> cfgRefR{"cfgRefR", 1.2, "reference R (in m) for extrapolation"}; // https://cds.cern.ch/record/1419204

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DielectronCut::PIDSchemes::kTOFif), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3, kTOFif : 4, kPIDML : 5]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    // Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    // Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -0.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +0.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -0.0, "min. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +0.0, "max. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -0.0, "min. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +0.0, "max. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_min_pin_pirejTPC{"cfg_min_pin_pirejTPC", 0.5, "min. pin for pion rejection in TPC"};
    Configurable<float> cfg_max_pin_pirejTPC{"cfg_max_pin_pirejTPC", 0.5, "max. pin for pion rejection in TPC"};
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

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};

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

  ~vpPairQC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms(&fRegistry);

    const AxisSpec axis_pt{ConfPtlBins, "p_{T,e} (GeV/c)"};
    const AxisSpec axis_eta{20, -1.0, +1.0, "#eta_{e}"};
    const AxisSpec axis_phi{36, 0.0, 2 * M_PI, "#varphi_{e} (rad.)"};
    const AxisSpec axis_dca{{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCA_{e}^{3D} (#sigma)"};
    fRegistry.add("Track/positive/hs", "rec. single electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_dca}, true);
    fRegistry.add("Track/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/positive/hPResolution", "p resolution;p (GeV/c);#Deltap/p", kTH2F, {{1000, 0.0f, 10.0f}, {100, 0.0f, 0.1f}}, false);
    fRegistry.add("Track/positive/hPtResolution", "p_{T} resolution;p (GeV/c);#Deltap_{T}/p_{T}", kTH2F, {{1000, 0.0f, 10.0f}, {100, 0.0f, 0.1f}}, false);
    fRegistry.add("Track/positive/hThetaResolution", "#theta resolution;p (GeV/c);#Delta#theta (rad.)", kTH2F, {{1000, 0.0f, 10.0f}, {100, 0.0f, 0.01f}}, false);
    fRegistry.add("Track/positive/hEtaResolution", "#eta resolution;p (GeV/c);#Delta#eta", kTH2F, {{1000, 0.0f, 10.0f}, {100, 0.0f, 0.01f}}, false);
    fRegistry.add("Track/positive/hPhiResolution", "#varphi resolution;p (GeV/c);#Delta#varphi (rad.)", kTH2F, {{1000, 0.0f, 10.0f}, {100, 0.0f, 0.01f}}, false);
    fRegistry.add("Track/positive/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/positive/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/positive/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
    fRegistry.add("Track/positive/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
    fRegistry.add("Track/positive/hDeltaPin", "p_{in} vs. p_{pv};p_{pv} (GeV/c);(p_{in} - p_{pv})/p_{pv}", kTH2F, {{1000, 0, 10}, {200, -1, +1}}, false);
    fRegistry.add("Track/positive/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/positive/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/positive/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/positive/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/positive/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/positive/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    fRegistry.add("Track/positive/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/positive/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/positive/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/positive/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/positive/hTOFbeta", "TOF #beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
    fRegistry.add("Track/positive/hChi2TOF", "TOF Chi2;p_{pv} (GeV/c);chi2", kTH2F, {{1000, 0, 10}, {100, 0, 10}}, false);
    fRegistry.add("Track/positive/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda);", kTH2F, {{1000, 0.f, 10.f}, {160, 0, 16}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    // fRegistry.add("Track/positive/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/positive/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    // fRegistry.add("Track/positive/hTOFNsigmaMu", "TOF n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    // fRegistry.add("Track/positive/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    // fRegistry.add("Track/positive/hTOFNsigmaKa", "TOF n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    // fRegistry.add("Track/positive/hTOFNsigmaPr", "TOF n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.addClone("Track/positive/", "Track/negative/");

    const AxisSpec axis_mass{50, 0, 0.05, "m_{ee} (GeV/c^{2})"};
    const AxisSpec axis_pair_pt{200, 0, 2, "p_{T,ee} (GeV/c)"};
    const AxisSpec axis_pair_dca_3d{100, 0, 10, "DCA_{ee}^{3D} (#sigma)"};
    const AxisSpec axis_pair_dca_xy{100, 0, 10, "DCA_{ee}^{XY} (#sigma)"};
    const AxisSpec axis_phiv{90, 0, M_PI, "#varphi_{V} (rad.)"};
    // for pair
    fRegistry.add("Pair/hMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {axis_mass, axis_pair_pt}, true);
    fRegistry.add("Pair/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {axis_phiv, axis_mass}, true);
    fRegistry.add("Pair/hDCA3DvsPhiV", "DCA_{ee}^{3D} vs. #varphi_{V};#varphi_{V} (rad.);DCA_{ee}^{3D} (#sigma)", kTH2F, {axis_phiv, axis_pair_dca_3d}, true);
    fRegistry.add("Pair/hDCAXYvsPhiV", "DCA_{ee}^{XY} vs. #varphi_{V};#varphi_{V} (rad.);DCA_{ee}^{XY} (#sigma)", kTH2F, {axis_phiv, axis_pair_dca_xy}, true);
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
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
  }

  o2::analysis::MlResponseDielectronSingleTrack<float> mlResponseSingleTrack;
  void DefineDielectronCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for pair
    fDielectronCut.SetMeeRange(dielectroncuts.cfg_min_mass, dielectroncuts.cfg_max_mass);
    fDielectronCut.SetPairPtRange(dielectroncuts.cfg_min_pair_pt, dielectroncuts.cfg_max_pair_pt);
    fDielectronCut.SetPairYRange(dielectroncuts.cfg_min_pair_y, dielectroncuts.cfg_max_pair_y);
    fDielectronCut.SetPairDCARange(dielectroncuts.cfg_min_pair_dca3d, dielectroncuts.cfg_max_pair_dca3d); // in sigma
    fDielectronCut.SetMaxMeePhiVDep([&](float phiv) { return dielectroncuts.cfg_phiv_intercept + phiv * dielectroncuts.cfg_phiv_slope; }, dielectroncuts.cfg_min_phiv, dielectroncuts.cfg_max_phiv);
    fDielectronCut.ApplyPhiV(dielectroncuts.cfg_apply_phiv);
    fDielectronCut.ApplyPrefilter(dielectroncuts.cfg_apply_pf);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, dielectroncuts.cfg_max_pt_track);
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetTrackDca3DRange(0.f, dielectroncuts.cfg_max_dca3dsigma_track); // in sigma
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
    fDielectronCut.SetChi2TOF(0.0, dielectroncuts.cfg_max_chi2tof);
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
    fDielectronCut.SetPinRangeForPionRejectionTPC(dielectroncuts.cfg_min_pin_pirejTPC, dielectroncuts.cfg_max_pin_pirejTPC);

    if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      static constexpr int nClassesMl = 2;
      const std::vector<int> cutDirMl = {o2::cuts_ml::CutSmaller, o2::cuts_ml::CutNot};
      const std::vector<std::string> labelsClasses = {"Signal", "Background"};
      const uint32_t nBinsMl = dielectroncuts.binsMl.value.size() - 1;
      const std::vector<std::string> labelsBins(nBinsMl, "bin");
      double cutsMlArr[nBinsMl][nClassesMl];
      for (uint32_t i = 0; i < nBinsMl; i++) {
        cutsMlArr[i][0] = dielectroncuts.cutsMl.value[i];
        cutsMlArr[i][1] = 0.;
      }
      o2::framework::LabeledArray<double> cutsMl = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      mlResponseSingleTrack.configure(dielectroncuts.binsMl.value, cutsMl, cutDirMl, nClassesMl);
      if (dielectroncuts.loadModelsFromCCDB) {
        ccdbApi.init(ccdburl);
        mlResponseSingleTrack.setModelPathsCCDB(dielectroncuts.onnxFileNames.value, ccdbApi, dielectroncuts.onnxPathsCCDB.value, dielectroncuts.timestampCCDB.value);
      } else {
        mlResponseSingleTrack.setModelPathsLocal(dielectroncuts.onnxFileNames.value);
      }
      mlResponseSingleTrack.cacheInputFeaturesIndices(dielectroncuts.namesInputFeatures);
      mlResponseSingleTrack.cacheBinningIndex(dielectroncuts.nameBinningFeature);
      mlResponseSingleTrack.init(dielectroncuts.enableOptimizations.value);

      fDielectronCut.SetPIDMlResponse(&mlResponseSingleTrack);
    } // end of PID ML
  }

  template <typename TCollision, typename TTrack1, typename TTrack2>
  bool fillPairInfo(TCollision const&, TTrack1 const& t1, TTrack2 const& t2)
  {
    // if (t1.trackId() == t2.trackId()) { // this is protection against pairing identical 2 tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
    //   return false;
    // }

    if (!fDielectronCut.IsSelectedTrack(t1) || !fDielectronCut.IsSelectedTrack(t2)) {
      return false;
    }

    if (!fDielectronCut.IsSelectedPair(t1, t2, d_bz, dielectroncuts.cfgRefR)) {
      return false;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);

    float dca_t1 = dca3DinSigma(t1);
    float dca_t2 = dca3DinSigma(t2);
    float pair_dca_3d = std::sqrt((dca_t1 * dca_t1 + dca_t2 * dca_t2) / 2.);
    float pair_dca_xy = std::sqrt((std::pow(t1.dcaXY() / std::sqrt(t1.cYY()), 2) + std::pow(t2.dcaXY() / std::sqrt(t2.cYY()), 2)) / 2.);

    fRegistry.fill(HIST("Pair/hMvsPt"), v12.M(), v12.Pt());
    fRegistry.fill(HIST("Pair/hMvsPhiV"), phiv, v12.M());
    fRegistry.fill(HIST("Pair/hDCA3DvsPhiV"), phiv, pair_dca_3d);
    fRegistry.fill(HIST("Pair/hDCAXYvsPhiV"), phiv, pair_dca_xy);

    if (std::find(used_trackIds.begin(), used_trackIds.end(), t1.globalIndex()) == used_trackIds.end()) {
      used_trackIds.emplace_back(t1.globalIndex());
      fillTrackInfo(t1);
    }
    if (std::find(used_trackIds.begin(), used_trackIds.end(), t2.globalIndex()) == used_trackIds.end()) {
      used_trackIds.emplace_back(t2.globalIndex());
      fillTrackInfo(t2);
    }

    return true;
  }

  template <typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    float weight = 1.f;
    float dca_3d = dca3DinSigma(track);
    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/positive/hs"), track.pt(), track.eta(), track.phi(), dca_3d, weight);
      fRegistry.fill(HIST("Track/positive/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/positive/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/positive/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/positive/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/positive/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/positive/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/positive/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/positive/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/positive/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/positive/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/positive/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/positive/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
      fRegistry.fill(HIST("Track/positive/hDeltaPin"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
      fRegistry.fill(HIST("Track/positive/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/positive/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/positive/hChi2TOF"), track.p(), track.tofChi2());

      fRegistry.fill(HIST("Track/positive/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/positive/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/positive/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      // fRegistry.fill(HIST("Track/positive/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      // fRegistry.fill(HIST("Track/positive/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
    } else {
      fRegistry.fill(HIST("Track/negative/hs"), track.pt(), track.eta(), track.phi(), dca_3d, weight);
      fRegistry.fill(HIST("Track/negative/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/negative/hDCAxyz"), track.dcaXY(), track.dcaZ());
      fRegistry.fill(HIST("Track/negative/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
      fRegistry.fill(HIST("Track/negative/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/negative/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/negative/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/negative/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/negative/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/negative/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/negative/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/negative/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/negative/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
      fRegistry.fill(HIST("Track/negative/hDeltaPin"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
      fRegistry.fill(HIST("Track/negative/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/negative/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/negative/hChi2TOF"), track.p(), track.tofChi2());

      fRegistry.fill(HIST("Track/negative/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/negative/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/negative/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      // fRegistry.fill(HIST("Track/negative/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      // fRegistry.fill(HIST("Track/negative/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
    }
  }

  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin < o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  SliceCache cache;
  Preslice<MyTracks> perCollision_track = aod::emprimaryelectron::emeventId;
  Filter trackFilter_electron = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter_electron = dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl;
  Filter ttcaFilter_electron = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  using FilteredMyTracks = soa::Filtered<MyTracks>;

  Partition<FilteredMyTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  std::vector<int> used_trackIds;
  void processQC(FilteredMyCollisions const& collisions, FilteredMyTracks const& /*tracks*/)
  {
    for (auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillPairInfo(collision, pos, ele);
      }
    } // end of collision loop

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();

  } // end of process
  PROCESS_SWITCH(vpPairQC, processQC, "run vp pair QC", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(vpPairQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<vpPairQC>(cfgc, TaskName{"vp-pair-qc"})};
}
