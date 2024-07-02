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
// Analysis task for dimuon.
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

#include "DCAFitter/FwdDCAFitterN.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/DimuonCut.h"
#include "PWGEM/PhotonMeson/Core/EMEventCut.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EMFwdTrack.h"
#include "PWGEM/PhotonMeson/Utils/EventMixingHandler.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::photonmeson::utils::emfwdtrack;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds>;
using MyTrack = MyTracks::iterator;

using MyEMH = o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int, int>, std::pair<int, int>, EMFwdTrackWithCov>;

struct dimuonQC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgQvecEstimator{"cfgQvecEstimator", 3, "FT0M:0, FT0A:1, FT0C:2, BTOT:3"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<bool> cfgDo_v2{"cfgDo_v2", false, "flag to analyze v2"};
  Configurable<bool> cfgDo_v3{"cfgDo_v3", false, "flag to analyze v3"};
  Configurable<float> minY{"minY", -4.0, "minimum rapidity for reconstructed pairs"};
  Configurable<float> maxY{"maxY", -2.5, "maximum rapidity for reconstructed pairs"};
  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {VARIABLE_WIDTH, -M_PI / 2, -M_PI / 4, 0.0f, +M_PI / 4, +M_PI / 2}, "Mixing bins - event plane angle"};
  ConfigurableAxis ConfOccupancyBins{"ConfOccupancyBins", {VARIABLE_WIDTH, -1, 1e+10}, "Mixing bins - occupancy"};

  ConfigurableAxis ConfMmumuBins{"ConfMmumuBins", {VARIABLE_WIDTH, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70, 5.80, 5.90, 6.00, 6.10, 6.20, 6.30, 6.40, 6.50, 6.60, 6.70, 6.80, 6.90, 7.00, 7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.80, 7.90, 8.00, 8.10, 8.20, 8.30, 8.40, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.00, 11.10, 11.20, 11.30, 11.40, 11.50, 11.60, 11.70, 11.80, 11.90, 12.00}, "mmumu bins for output histograms"};
  ConfigurableAxis ConfPtmumuBins{"ConfPtmumuBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00}, "pTmumu bins for output histograms"};
  ConfigurableAxis ConfDCAmumuBins{"ConfDCAmumuBins", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAmumu bins for output histograms"};
  // ConfigurableAxis ConfPCAmumuBins{"ConfPCAmumuBins", {VARIABLE_WIDTH, 0.0, 0.5, 1, 1.5, 2, 3, 4, 5}, "PCAmumu bins for output histograms in mm"};

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

  DimuonCut fDimuonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dimuoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pt"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pt"};
    Configurable<float> cfg_min_pair_dcaxy{"cfg_min_pair_dcaxy", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dcaxy{"cfg_max_pair_dcaxy", 1e+10, "max pair dca3d in sigma"};

    Configurable<int> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
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
  } dimuoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;
  o2::vertexing::FwdDCAFitterN<2> fitter;
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

    emh_pos = new MyEMH(ndepth);
    emh_neg = new MyEMH(ndepth);

    DefineEMEventCut();
    DefineDimuonCut();
    addhistograms();

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(90.f);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setTGeoMat(false);
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

    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
  }

  ~dimuonQC()
  {
    delete emh_pos;
    emh_pos = 0x0;
    delete emh_neg;
    emh_neg = 0x0;

    map_mixed_eventId_to_q2vector.clear();
    map_mixed_eventId_to_q3vector.clear();

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();
  }

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry, cfgDo_v2 | cfgDo_v3);

    // pair info
    const AxisSpec axis_mass{ConfMmumuBins, "m_{#mu#mu} (GeV/c^{2})"};
    const AxisSpec axis_pt{ConfPtmumuBins, "p_{T,#mu#mu} (GeV/c)"};
    const AxisSpec axis_dca{ConfDCAmumuBins, "DCA_{#mu#mu}^{xy} (#sigma)"};
    // const AxisSpec axis_pca{ConfPCAmumuBins, "PCA (mm)"}; // particle closest approach

    std::string_view qvec_det_names[4] = {"FT0M", "FT0A", "FT0C", "BTOT"};
    int nbin_sp2 = 1;
    int nbin_sp3 = 1;
    if (cfgDo_v2) {
      nbin_sp2 = 100;
    }
    if (cfgDo_v3) {
      nbin_sp3 = 100;
    }
    const AxisSpec axis_sp2{nbin_sp2, -5.f, 5.f, Form("u_{2}^{#mu#mu} #upoint Q_{2}^{%s}", qvec_det_names[cfgQvecEstimator].data())};
    const AxisSpec axis_sp3{nbin_sp3, -5.f, 5.f, Form("u_{3}^{#mu#mu} #upoint Q_{3}^{%s}", qvec_det_names[cfgQvecEstimator].data())};

    fRegistry.add("Pair/same/uls/hs", "dimuon", kTHnSparseD, {axis_mass, axis_pt, axis_dca, axis_sp2, axis_sp3}, true);
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");
    fRegistry.addClone("Pair/same/", "Pair/mix/");

    // for track info
    fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {25, -4.5f, -2.0f}}, false);
    fRegistry.add("Track/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
    fRegistry.add("Track/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
    fRegistry.add("Track/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
    fRegistry.add("Track/hPDCA", "pDCA;p_{T} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
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

  template <int ev_id, typename TCollision, typename TTrack1, typename TTrack2>
  bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, const float q2x_mixed = 0.f, const float q2y_mixed = 0.f, const float q3x_mixed = 0.f, const float q3y_mixed = 0.f)
  {
    if constexpr (ev_id == 1) {
      if (t1.has_ambiguousMuons() && t2.has_ambiguousMuons()) {
        for (auto& possible_id1 : t1.ambiguousMuonsIds()) {
          for (auto& possible_id2 : t2.ambiguousMuonsIds()) {
            if (possible_id1 == possible_id2) {
              // LOGF(info, "event id = %d: same track is found. t1.fwdtrackId() = %d, t1.collisionId() = %d, t1.pt() = %f, t1.eta() = %f, t1.phi() = %f, t2.fwdtrackId() = %d, t2.collisionId() = %d, t2.pt() = %f, t2.eta() = %f, t2.phi() = %f", ev_id, t1.fwdtrackId(), t1.collisionId(), t1.pt(), t1.eta(), t1.phi(), t2.fwdtrackId(), t2.collisionId(), t2.pt(), t2.eta(), t2.phi());
              return false; // this is protection against pairing 2 identical tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
            }
          }
        }
      }
    }

    if constexpr (ev_id == 0) {
      if (!fDimuonCut.IsSelectedTrack(t1) || !fDimuonCut.IsSelectedTrack(t2)) {
        return false;
      }
    }

    if (!fDimuonCut.IsSelectedPair(t1, t2)) {
      return false;
    }

    // float pca = 999.f, lxy = 999.f; // in unit of cm
    // o2::aod::pwgem::dilepton::utils::pairutil::isSVFoundFwd(fitter, collision, t1, t2, pca, lxy);

    std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
    std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
    std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
    std::array<float, 2> q2btot = {collision.q2xbtot(), collision.q2ybtot()};
    const std::array<float, 2> q2vector[4] = {q2ft0m, q2ft0a, q2ft0c, q2btot};

    std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
    std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
    std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
    std::array<float, 2> q3btot = {collision.q3xbtot(), collision.q3ybtot()};
    const std::array<float, 2> q3vector[4] = {q3ft0m, q3ft0a, q3ft0c, q3btot};

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    if (v12.Rapidity() < minY || maxY < v12.Rapidity()) {
      return false;
    }

    float dca_xy_t1 = fwdDcaXYinSigma(t1);
    float dca_xy_t2 = fwdDcaXYinSigma(t2);
    float dca_mumu_xy = std::sqrt((dca_xy_t1 * dca_xy_t1 + dca_xy_t2 * dca_xy_t2) / 2.);

    float sp2 = 0.f, sp3 = 0.f;
    if constexpr (ev_id == 0) {
      if (cfgDo_v2) {
        std::array<float, 2> u2_ll = {static_cast<float>(std::cos(2 * v12.Phi())), static_cast<float>(std::sin(2 * v12.Phi()))};
        sp2 = RecoDecay::dotProd(u2_ll, q2vector[cfgQvecEstimator]);
      }
      if (cfgDo_v3) {
        std::array<float, 2> u3_ll = {static_cast<float>(std::cos(3 * v12.Phi())), static_cast<float>(std::sin(3 * v12.Phi()))};
        sp3 = RecoDecay::dotProd(u3_ll, q3vector[cfgQvecEstimator]);
      }
    } else if constexpr (ev_id == 1) {
      if (cfgDo_v2) {
        std::array<float, 2> u2_l1 = {static_cast<float>(std::cos(2 * v1.Phi())), static_cast<float>(std::sin(2 * v1.Phi()))};
        std::array<float, 2> u2_l2 = {static_cast<float>(std::cos(2 * v2.Phi())), static_cast<float>(std::sin(2 * v2.Phi()))};
        sp2 = RecoDecay::dotProd(u2_l1, std::array<float, 2>{q2vector[cfgQvecEstimator][0], q2vector[cfgQvecEstimator][1]}) * std::cos(2 * (v1.Phi() - v12.Phi())) + RecoDecay::dotProd(u2_l2, std::array<float, 2>{q2x_mixed, q2y_mixed}) * std::cos(2 * (v2.Phi() - v12.Phi()));
      }
      if (cfgDo_v3) {
        std::array<float, 2> u3_l1 = {static_cast<float>(std::cos(3 * v1.Phi())), static_cast<float>(std::sin(3 * v1.Phi()))};
        std::array<float, 2> u3_l2 = {static_cast<float>(std::cos(3 * v2.Phi())), static_cast<float>(std::sin(3 * v2.Phi()))};
        sp3 = RecoDecay::dotProd(u3_l1, std::array<float, 2>{q3vector[cfgQvecEstimator][0], q3vector[cfgQvecEstimator][1]}) * std::cos(3 * (v1.Phi() - v12.Phi())) + RecoDecay::dotProd(u3_l2, std::array<float, 2>{q3x_mixed, q3y_mixed}) * std::cos(3 * (v2.Phi() - v12.Phi()));
      }
    }

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("uls/hs"), v12.M(), v12.Pt(), dca_mumu_xy, sp2, sp3);
    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lspp/hs"), v12.M(), v12.Pt(), dca_mumu_xy, sp2, sp3);
    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("lsmm/hs"), v12.M(), v12.Pt(), dca_mumu_xy, sp2, sp3);
    }

    // store tracks for event mixing without double counting
    if constexpr (ev_id == 0) {
      std::pair<int, int> key_df_collision = std::make_pair(ndf, collision.globalIndex());
      std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, t1.globalIndex());
      std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, t2.globalIndex());

      std::vector<int> possibleIds1;
      std::vector<int> possibleIds2;
      std::copy(t1.ambiguousMuonsIds().begin(), t1.ambiguousMuonsIds().end(), std::back_inserter(possibleIds1));
      std::copy(t2.ambiguousMuonsIds().begin(), t2.ambiguousMuonsIds().end(), std::back_inserter(possibleIds2));

      if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id1) == used_trackIds.end()) {
        used_trackIds.emplace_back(pair_tmp_id1);
        fillTrackInfo(t1);
        if (cfgDoMix) {
          if (t1.sign() > 0) {
            emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(t1.globalIndex(), collision.globalIndex(), t1.fwdtrackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), possibleIds1,
                                                                             t1.x(), t1.y(), t1.z(), t1.tgl(), t1.cXX(), t1.cXY(), t1.cYY(),
                                                                             t1.cPhiX(), t1.cPhiY(), t1.cPhiPhi(), t1.cTglX(), t1.cTglY(), t1.cTglPhi(), t1.cTglTgl(), t1.c1PtX(), t1.c1PtY(), t1.c1PtPhi(), t1.c1PtTgl(), t1.c1Pt21Pt2(), t1.chi2()));
          } else {
            emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(t1.globalIndex(), collision.globalIndex(), t1.fwdtrackId(), t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon, t1.sign(), t1.fwdDcaX(), t1.fwdDcaY(), possibleIds1,
                                                                             t1.x(), t1.y(), t1.z(), t1.tgl(), t1.cXX(), t1.cXY(), t1.cYY(),
                                                                             t1.cPhiX(), t1.cPhiY(), t1.cPhiPhi(), t1.cTglX(), t1.cTglY(), t1.cTglPhi(), t1.cTglTgl(), t1.c1PtX(), t1.c1PtY(), t1.c1PtPhi(), t1.c1PtTgl(), t1.c1Pt21Pt2(), t1.chi2()));
          }
        }
      }
      if (std::find(used_trackIds.begin(), used_trackIds.end(), pair_tmp_id2) == used_trackIds.end()) {
        used_trackIds.emplace_back(pair_tmp_id2);
        fillTrackInfo(t2);
        if (cfgDoMix) {
          if (t2.sign() > 0) {
            emh_pos->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                             t2.x(), t2.y(), t2.z(), t2.tgl(), t2.cXX(), t2.cXY(), t2.cYY(),
                                                                             t2.cPhiX(), t2.cPhiY(), t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(), t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2(), t2.chi2()));
          } else {
            emh_neg->AddTrackToEventPool(key_df_collision, EMFwdTrackWithCov(t2.globalIndex(), collision.globalIndex(), t2.fwdtrackId(), t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon, t2.sign(), t2.fwdDcaX(), t2.fwdDcaY(), possibleIds2,
                                                                             t2.x(), t2.y(), t2.z(), t2.tgl(), t2.cXX(), t2.cXY(), t2.cYY(),
                                                                             t2.cPhiX(), t2.cPhiY(), t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(), t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2(), t2.chi2()));
          }
        }
      }
    }
    return true;
  }

  template <typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    fRegistry.fill(HIST("Track/hPt"), track.pt());
    fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/hEtaPhi"), track.phi(), track.eta());
    fRegistry.fill(HIST("Track/hTrackType"), track.trackType());
    fRegistry.fill(HIST("Track/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
    fRegistry.fill(HIST("Track/hDCAxySigma"), track.fwdDcaX() / sqrt(track.cXX()), track.fwdDcaY() / sqrt(track.cYY()));
    fRegistry.fill(HIST("Track/hDCAxRes_Pt"), track.pt(), sqrt(track.cXX()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/hDCAyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/hNclsMFT"), track.nClustersMFT());
    fRegistry.fill(HIST("Track/hNclsMCH"), track.nClusters());
    fRegistry.fill(HIST("Track/hPDCA"), track.pt(), track.pDca());
    fRegistry.fill(HIST("Track/hChi2"), track.chi2());
    fRegistry.fill(HIST("Track/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
    fRegistry.fill(HIST("Track/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
    fRegistry.fill(HIST("Track/hMFTClusterMap"), track.mftClusterMap());
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  SliceCache cache;
  Preslice<MyTracks> perCollision_track = aod::emprimarymuon::emeventId;
  // Filter trackFilter =  o2::aod::fwdtrack::pt > 1.f * dimuoncuts.cfg_min_pt_track.node();
  Filter trackFilter = static_cast<float>(dimuoncuts.cfg_min_pt_track) < o2::aod::fwdtrack::pt && static_cast<float>(dimuoncuts.cfg_min_eta_track) < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < static_cast<float>(dimuoncuts.cfg_max_eta_track);
  // Filter trackFilter = o2::aod::fwdtrack::trackType == static_cast<uint8_t>(dimuoncuts.cfg_track_type) && static_cast<float>(dimuoncuts.cfg_min_pt_track) < o2::aod::fwdtrack::pt && static_cast<float>(dimuoncuts.cfg_min_eta_track) < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < static_cast<float>(dimuoncuts.cfg_max_eta_track);
  using FilteredMyTracks = soa::Filtered<MyTracks>;

  MyEMH* emh_pos = nullptr;
  MyEMH* emh_neg = nullptr;

  Partition<FilteredMyTracks> posTracks = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<FilteredMyTracks> negTracks = o2::aod::emprimarymuon::sign < int8_t(0);

  std::map<std::pair<int, int>, std::array<float, 2>> map_mixed_eventId_to_q2vector;
  std::map<std::pair<int, int>, std::array<float, 2>> map_mixed_eventId_to_q3vector;

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

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, cfgDo_v2 | cfgDo_v3);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, cfgDo_v2 | cfgDo_v3);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2btot = {collision.q2xbtot(), collision.q2ybtot()};
      const std::array<float, 2> q2vector[4] = {q2ft0m, q2ft0a, q2ft0c, q2btot};

      std::array<float, 2> q3ft0m = {collision.q3xft0m(), collision.q3yft0m()};
      std::array<float, 2> q3ft0a = {collision.q3xft0a(), collision.q3yft0a()};
      std::array<float, 2> q3ft0c = {collision.q3xft0c(), collision.q3yft0c()};
      std::array<float, 2> q3btot = {collision.q3xbtot(), collision.q3ybtot()};
      const std::array<float, 2> q3vector[4] = {q3ft0m, q3ft0a, q3ft0c, q3btot};

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      // LOGF(info, "collision.globalIndex() = %d , collision.posZ() = %f , collision.numContrib() = %d, centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", collision.globalIndex(), collision.posZ(), collision.numContrib(), centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      // for (auto& neg : negTracks_per_coll) {
      //   LOGF(info, "neg.trackType() = %d , neg.pt() = %f", neg.trackType(), neg.pt());
      // }

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

      float ep2 = collision.ep2ft0c();
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
      auto selected_negTracks_in_this_event = emh_neg->GetTracksPerCollision(key_df_collision);
      // LOGF(info, "N selected tracks in current event (%d, %d), zvtx = %f, centrality = %f , npos = %d , nele = %d, nuls = %d , nlspp = %d, nlsmm = %d", ndf, collision.globalIndex(), collision.posZ(), centralities[cfgCentEstimator], selected_posTracks_in_this_event.size(), selected_negTracks_in_this_event.size(), nuls, nlspp, nlsmm);

      auto collisionIds_in_mixing_pool = emh_pos->GetCollisionIdsFromEventPool(key_bin); // pos/ele does not matter.
      // LOGF(info, "collisionIds_in_mixing_pool.size() = %d", collisionIds_in_mixing_pool.size());

      for (auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
        int mix_dfId = mix_dfId_collisionId.first;
        int mix_collisionId = mix_dfId_collisionId.second;
        if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
          continue;
        }

        float q2x_mixed = map_mixed_eventId_to_q2vector[mix_dfId_collisionId][0];
        float q2y_mixed = map_mixed_eventId_to_q2vector[mix_dfId_collisionId][1];
        float q3x_mixed = map_mixed_eventId_to_q3vector[mix_dfId_collisionId][0];
        float q3y_mixed = map_mixed_eventId_to_q3vector[mix_dfId_collisionId][1];

        auto posTracks_from_event_pool = emh_pos->GetTracksPerCollision(mix_dfId_collisionId);
        auto negTracks_from_event_pool = emh_neg->GetTracksPerCollision(mix_dfId_collisionId);
        // LOGF(info, "Do event mixing: current event (%d, %d) | event pool (%d, %d), npos = %d , nele = %d", ndf, collision.globalIndex(), mix_dfId, mix_collisionId, posTracks_from_event_pool.size(), negTracks_from_event_pool.size());

        for (auto& pos : selected_posTracks_in_this_event) { // ULS mix
          for (auto& ele : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos, ele, q2x_mixed, q2y_mixed, q3x_mixed, q3y_mixed);
          }
        }

        for (auto& ele : selected_negTracks_in_this_event) { // ULS mix
          for (auto& pos : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, ele, pos, q2x_mixed, q2y_mixed, q3x_mixed, q3y_mixed);
          }
        }

        for (auto& pos1 : selected_posTracks_in_this_event) { // LS++ mix
          for (auto& pos2 : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos1, pos2, q2x_mixed, q2y_mixed, q3x_mixed, q3y_mixed);
          }
        }

        for (auto& ele1 : selected_negTracks_in_this_event) { // LS-- mix
          for (auto& ele2 : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, ele1, ele2, q2x_mixed, q2y_mixed, q3x_mixed, q3y_mixed);
          }
        }
      } // end of loop over mixed event pool

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) {
        map_mixed_eventId_to_q2vector[key_df_collision] = q2vector[cfgQvecEstimator];
        map_mixed_eventId_to_q3vector[key_df_collision] = q3vector[cfgQvecEstimator];
        emh_pos->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh_neg->AddCollisionIdAtLast(key_bin, key_df_collision);
      }

    } // end of collision loop

    ndf++;
  } // end of process
  PROCESS_SWITCH(dimuonQC, processQC, "run dimuon QC", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(dimuonQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dimuonQC>(cfgc, TaskName{"dimuon-qc"})};
}
