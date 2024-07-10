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

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "DCAFitter/FwdDCAFitterN.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
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

using MyMCTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMPrimaryMuonMCLabels>;
using MyMCTrack = MyMCTracks::iterator;

struct dimuonQCMC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> minY{"minY", -4.0, "minimum rapidity for reconstructed pairs"};
  Configurable<float> maxY{"maxY", -2.5, "maximum rapidity for reconstructed pairs"};

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
  o2::vertexing::FwdDCAFitterN<2> fitter;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int mRunNumber;
  float d_bz;

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt{"min_mcPt", 0.05, "min. MC pT"};
    Configurable<float> max_mcPt{"max_mcPt", 1e+10, "max. MC pT"};
    Configurable<float> min_mcEta{"min_mcEta", -4.0, "min. MC eta"};
    Configurable<float> max_mcEta{"max_mcEta", -2.5, "max. MC eta"};
  } mctrackcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view ele_source_types[8] = {"lf/", "PromptJPsi/", "NonPromptJPsi/", "PromptPsi2S/", "NonPromptPsi2S/", "c2mu/", "b2mu/", "b2c2mu/"};

  ~dimuonQCMC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms(&fRegistry, cfgDoFlow);

    const AxisSpec axis_mass{ConfMmumuBins, "m_{#mu#mu} (GeV/c^{2})"};
    const AxisSpec axis_pt{ConfPtmumuBins, "p_{T,#mu#mu} (GeV/c)"};
    const AxisSpec axis_dca{ConfDCAmumuBins, "DCA_{#mu#mu}^{xy} (#sigma)"};
    // const AxisSpec axis_pca{ConfPCAmumuBins, "PCA (mm)"}; // particle closest approach

    // generated info
    fRegistry.add("Generated/sm/Eta/hMvsPt", "m_{#mu#mu} vs. p_{T,#mu#mu} ULS", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/EtaPrime/");
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/Rho/");
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/Omega/");
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/Phi/");
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/PromptJPsi/");
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/NonPromptJPsi/");
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/PromptPsi2S/");
    fRegistry.addClone("Generated/sm/Eta/", "Generated/sm/NonPromptPsi2S/");

    fRegistry.add("Generated/ccbar/c2mu_c2mu/hadron_hadron/hMvsPt", "m_{#mu#mu} vs. p_{T,#mu#mu}", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry.addClone("Generated/ccbar/c2mu_c2mu/hadron_hadron/", "Generated/ccbar/c2mu_c2mu/meson_meson/");
    fRegistry.addClone("Generated/ccbar/c2mu_c2mu/hadron_hadron/", "Generated/ccbar/c2mu_c2mu/baryon_baryon/");
    fRegistry.addClone("Generated/ccbar/c2mu_c2mu/hadron_hadron/", "Generated/ccbar/c2mu_c2mu/meson_baryon/");
    fRegistry.addClone("Generated/ccbar/c2mu_c2mu/", "Generated/bbbar/b2mu_b2mu/");
    fRegistry.addClone("Generated/ccbar/c2mu_c2mu/", "Generated/bbbar/b2c2mu_b2c2mu/");
    fRegistry.addClone("Generated/ccbar/c2mu_c2mu/", "Generated/bbbar/b2c2mu_b2mu_sameb/");
    fRegistry.addClone("Generated/ccbar/c2mu_c2mu/", "Generated/bbbar/b2c2mu_b2mu_diffb/"); // LS

    // reconstructed pair info
    fRegistry.add("Pair/sm/Eta/hs", "hs pair", kTHnSparseF, {axis_mass, axis_pt, axis_dca}, true);
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/EtaPrime/");
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/Rho/");
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/Omega/");
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/Phi/");
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/PromptJPsi/");
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/NonPromptJPsi/");
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/PromptPsi2S/");
    fRegistry.addClone("Pair/sm/Eta/", "Pair/sm/NonPromptPsi2S/");

    fRegistry.add("Pair/ccbar/c2mu_c2mu/hadron_hadron/hs", "hs pair", kTHnSparseF, {axis_mass, axis_pt, axis_dca}, true);
    fRegistry.addClone("Pair/ccbar/c2mu_c2mu/hadron_hadron/", "Pair/ccbar/c2mu_c2mu/meson_meson/");
    fRegistry.addClone("Pair/ccbar/c2mu_c2mu/hadron_hadron/", "Pair/ccbar/c2mu_c2mu/baryon_baryon/");
    fRegistry.addClone("Pair/ccbar/c2mu_c2mu/hadron_hadron/", "Pair/ccbar/c2mu_c2mu/meson_baryon/");
    fRegistry.addClone("Pair/ccbar/c2mu_c2mu/", "Pair/bbbar/b2mu_b2mu/");
    fRegistry.addClone("Pair/ccbar/c2mu_c2mu/", "Pair/bbbar/b2c2mu_b2c2mu/");
    fRegistry.addClone("Pair/ccbar/c2mu_c2mu/", "Pair/bbbar/b2c2mu_b2mu_sameb/");
    fRegistry.addClone("Pair/ccbar/c2mu_c2mu/", "Pair/bbbar/b2c2mu_b2mu_diffb/"); // LS

    // track info
    fRegistry.add("Track/lf/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/lf/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/lf/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {30, -4.5f, -2.0f}}, false);
    fRegistry.add("Track/lf/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
    fRegistry.add("Track/lf/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/lf/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/lf/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/lf/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/lf/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
    fRegistry.add("Track/lf/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
    fRegistry.add("Track/lf/hPDCA", "pDCA;p_{T} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/lf/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/lf/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/lf/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/lf/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
    fRegistry.add("Track/lf/hPtGen_DeltaPtOverPtGen", "muon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/lf/hPtGen_DeltaEta", "muon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/lf/hPtGen_DeltaPhi", "muon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.addClone("Track/lf/", "Track/PromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/PromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/c2mu/");
    fRegistry.addClone("Track/lf/", "Track/b2mu/");
    fRegistry.addClone("Track/lf/", "Track/b2c2mu/");
  }

  bool cfgDoFlow = false;
  void init(InitContext&)
  {
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

  template <typename TTrack, typename TMCParticles>
  int FindLF(TTrack const& posmc, TTrack const& elemc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 221, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 331, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 113, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 223, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 333, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 100443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 100553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 200553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 300553, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
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
    // fill track info that belong to true pairs.
    if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
      auto mctrack = track.template emmcparticle_as<TMCParticles>();
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPt"), track.pt());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hEtaPhi"), track.phi(), track.eta());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hTrackType"), track.trackType());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXX()), track.fwdDcaY() / std::sqrt(track.cYY()));
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXX()) * 1e+4);
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4);
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNclsMCH"), track.nClusters());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hNclsMFT"), track.nClustersMFT());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPDCA"), track.pt(), track.pDca());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hChi2"), track.chi2());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hChi2MatchMCHMID"), track.chi2MatchMCHMID());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hMFTClusterMap"), track.mftClusterMap());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
      fRegistry.fill(HIST("Track/") + HIST(ele_source_types[e_source_id]) + HIST("hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());
      used_trackIds.emplace_back(track.globalIndex());
    }
  }

  template <typename TCollision, typename TTrack1, typename TTrack2, typename TMCParticles>
  bool fillTruePairInfo(TCollision const& /*collision*/, TTrack1 const& t1, TTrack2 const& t2, TMCParticles const& mcparticles)
  {
    if (!fDimuonCut.IsSelectedTrack(t1) || !fDimuonCut.IsSelectedTrack(t2)) {
      return false;
    }

    if (!fDimuonCut.IsSelectedPair(t1, t2)) {
      return false;
    }

    // float pca = 999.f, lxy = 999.f; // in unit of cm
    // o2::aod::pwgem::dilepton::utils::pairutil::isSVFoundFwd(fitter, collision, t1, t2, pca, lxy);

    auto t1mc = t1.template emmcparticle_as<TMCParticles>();
    auto t2mc = t2.template emmcparticle_as<TMCParticles>();

    if (abs(t1mc.pdgCode()) != 13 || abs(t2mc.pdgCode()) != 13) {
      return false;
    }

    int mother_id = FindLF(t1mc, t2mc, mcparticles);
    int hfee_type = IsHF(t1mc, t2mc, mcparticles);
    if (mother_id < 0 && hfee_type < 0) {
      return false;
    }
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    if (v12.Rapidity() < minY || maxY < v12.Rapidity()) {
      return false;
    }

    float dca_xy_t1 = fwdDcaXYinSigma(t1);
    float dca_xy_t2 = fwdDcaXYinSigma(t2);
    float dca_mumu_xy = std::sqrt((dca_xy_t1 * dca_xy_t1 + dca_xy_t2 * dca_xy_t2) / 2.);

    if (mother_id > -1 && t1mc.pdgCode() * t2mc.pdgCode() < 0) {
      auto mcmother = mcparticles.iteratorAt(mother_id);
      if (mcmother.isPhysicalPrimary() || mcmother.producedByGenerator()) {
        if ((t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && (t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
          switch (abs(mcmother.pdgCode())) {
            case 221:
              fRegistry.fill(HIST("Pair/sm/Eta/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              fillTrackInfo<0, TMCParticles>(t1);
              fillTrackInfo<0, TMCParticles>(t2);
              break;
            case 331:
              fRegistry.fill(HIST("Pair/sm/EtaPrime/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              fillTrackInfo<0, TMCParticles>(t1);
              fillTrackInfo<0, TMCParticles>(t2);
              break;
            case 113:
              fRegistry.fill(HIST("Pair/sm/Rho/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              fillTrackInfo<0, TMCParticles>(t1);
              fillTrackInfo<0, TMCParticles>(t2);
              break;
            case 223:
              fRegistry.fill(HIST("Pair/sm/Omega/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              fillTrackInfo<0, TMCParticles>(t1);
              fillTrackInfo<0, TMCParticles>(t2);
              break;
            case 333:
              fRegistry.fill(HIST("Pair/sm/Phi/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              fillTrackInfo<0, TMCParticles>(t1);
              fillTrackInfo<0, TMCParticles>(t2);
              break;
            case 443: {
              if (IsFromBeauty(mcmother, mcparticles) > 0) {
                fRegistry.fill(HIST("Pair/sm/NonPromptJPsi/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<2, TMCParticles>(t1);
                fillTrackInfo<2, TMCParticles>(t2);
              } else {
                fRegistry.fill(HIST("Pair/sm/PromptJPsi/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<1, TMCParticles>(t1);
                fillTrackInfo<1, TMCParticles>(t2);
              }
              break;
            }
            case 100443: {
              if (IsFromBeauty(mcmother, mcparticles) > 0) {
                fRegistry.fill(HIST("Pair/sm/NonPromptPsi2S/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<4, TMCParticles>(t1);
                fillTrackInfo<4, TMCParticles>(t2);
              } else {
                fRegistry.fill(HIST("Pair/sm/PromptPsi2S/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<3, TMCParticles>(t1);
                fillTrackInfo<3, TMCParticles>(t2);
              }
              break;
            }
            default:
              break;
          }
        }
      } // end of primary selection for same mother
    } else if (hfee_type > -1) {
      if ((t1mc.isPhysicalPrimary() || t1mc.producedByGenerator()) && (t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
        auto mp1 = mcparticles.iteratorAt(t1mc.mothersIds()[0]);
        auto mp2 = mcparticles.iteratorAt(t2mc.mothersIds()[0]);
        if (t1mc.pdgCode() * t2mc.pdgCode() < 0) { // ULS
          switch (hfee_type) {
            case static_cast<int>(EM_HFeeType::kCe_Ce): {
              fRegistry.fill(HIST("Pair/ccbar/c2mu_c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Pair/ccbar/c2mu_c2mu/meson_meson/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<5, TMCParticles>(t1);
                fillTrackInfo<5, TMCParticles>(t2);
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Pair/ccbar/c2mu_c2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<5, TMCParticles>(t1);
                fillTrackInfo<5, TMCParticles>(t2);
              } else {
                fRegistry.fill(HIST("Pair/ccbar/c2mu_c2mu/meson_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<5, TMCParticles>(t1);
                fillTrackInfo<5, TMCParticles>(t2);
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBe_Be): {
              fRegistry.fill(HIST("Pair/bbbar/b2mu_b2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              if (isBeautyMeson(mp1) && isBeautyMeson(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2mu_b2mu/meson_meson/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<6, TMCParticles>(t1);
                fillTrackInfo<6, TMCParticles>(t2);
              } else if (isBeautyBaryon(mp1) && isBeautyBaryon(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2mu_b2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<6, TMCParticles>(t1);
                fillTrackInfo<6, TMCParticles>(t2);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2mu_b2mu/meson_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<6, TMCParticles>(t1);
                fillTrackInfo<6, TMCParticles>(t2);
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_BCe): {
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2c2mu/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2c2mu/meson_meson/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<7, TMCParticles>(t1);
                fillTrackInfo<7, TMCParticles>(t2);
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2c2mu/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<7, TMCParticles>(t1);
                fillTrackInfo<7, TMCParticles>(t2);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2c2mu/meson_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
                fillTrackInfo<7, TMCParticles>(t1);
                fillTrackInfo<7, TMCParticles>(t2);
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): { // ULS
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_sameb/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_sameb/meson_meson/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_sameb/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_sameb/meson_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              }
              if ((isCharmMeson(mp1) || isCharmBaryon(mp1)) && (isBeautyMeson(mp2) || isBeautyBaryon(mp2))) {
                fillTrackInfo<6, TMCParticles>(t1);
                fillTrackInfo<7, TMCParticles>(t2);
              } else {
                fillTrackInfo<7, TMCParticles>(t1);
                fillTrackInfo<6, TMCParticles>(t2);
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
              fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_diffb/hadron_hadron/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_diffb/meson_meson/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_diffb/baryon_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              } else {
                fRegistry.fill(HIST("Pair/bbbar/b2c2mu_b2mu_diffb/meson_baryon/hs"), v12.M(), v12.Pt(), dca_mumu_xy);
              }
              if ((isCharmMeson(mp1) || isCharmBaryon(mp1)) && (isBeautyMeson(mp2) || isBeautyBaryon(mp2))) {
                fillTrackInfo<6, TMCParticles>(t1);
                fillTrackInfo<7, TMCParticles>(t2);
              } else {
                fillTrackInfo<7, TMCParticles>(t1);
                fillTrackInfo<6, TMCParticles>(t2);
              }
              break;
            }
            default:
              break;
          }
        }
      }
    } // end of HF evaluation
    return true;
  }

  std::vector<int> used_trackIds;
  SliceCache cache;
  Preslice<MyMCTracks> perCollision_track = aod::emprimarymuon::emeventId;
  Filter trackFilter = static_cast<float>(dimuoncuts.cfg_min_pt_track) < o2::aod::fwdtrack::pt && static_cast<float>(dimuoncuts.cfg_min_eta_track) < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < static_cast<float>(dimuoncuts.cfg_max_eta_track);

  using FilteredMyMCTracks = soa::Filtered<MyMCTracks>;
  Partition<FilteredMyMCTracks> posTracks = o2::aod::emprimarymuon::sign > int8_t(0);
  Partition<FilteredMyMCTracks> negTracks = o2::aod::emprimarymuon::sign < int8_t(0);

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

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, cfgDoFlow);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, cfgDoFlow);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillTruePairInfo(collision, pos, ele, mcparticles);
      } // end of ULS pair loop

      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        fillTruePairInfo(collision, pos1, pos2, mcparticles);
      } // end of ULS pair loop

      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS__
        fillTruePairInfo(collision, ele1, ele2, mcparticles);
      } // end of ULS pair loop

    } // end of collision loop

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();

  } // end of process
  PROCESS_SWITCH(dimuonQCMC, processQCMC, "run dimuon QC MC", true);

  Partition<aod::EMMCParticles> posTracksMC = o2::aod::mcparticle::pdgCode == -13; // mu+
  Partition<aod::EMMCParticles> negTracksMC = o2::aod::mcparticle::pdgCode == +13; // mu-
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

      auto posTracks_per_coll = posTracksMC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);

      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isInAcceptance(t1) || !isInAcceptance(t2)) {
          continue;
        }

        if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
          continue;
        }
        if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
          continue;
        }

        int mother_id = FindLF(t1, t2, mcparticles);
        int hfee_type = IsHF(t1, t2, mcparticles);
        if (mother_id < 0 && hfee_type < 0) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        if (v12.Rapidity() < minY || maxY < v12.Rapidity()) {
          continue;
        }

        if (mother_id > -1) {
          auto mcmother = mcparticles.iteratorAt(mother_id);
          if (mcmother.isPhysicalPrimary() || mcmother.producedByGenerator()) {

            switch (abs(mcmother.pdgCode())) {
              case 221:
                fRegistry.fill(HIST("Generated/sm/Eta/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 331:
                fRegistry.fill(HIST("Generated/sm/EtaPrime/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 113:
                fRegistry.fill(HIST("Generated/sm/Rho/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 223:
                fRegistry.fill(HIST("Generated/sm/Omega/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 333:
                fRegistry.fill(HIST("Generated/sm/Phi/hMvsPt"), v12.M(), v12.Pt());
                break;
              case 443: {
                if (IsFromBeauty(mcmother, mcparticles) > 0) {
                  fRegistry.fill(HIST("Generated/sm/NonPromptJPsi/hMvsPt"), v12.M(), v12.Pt());
                } else {
                  fRegistry.fill(HIST("Generated/sm/PromptJPsi/hMvsPt"), v12.M(), v12.Pt());
                }
                break;
              }
              case 100443: {
                if (IsFromBeauty(mcmother, mcparticles) > 0) {
                  fRegistry.fill(HIST("Generated/sm/NonPromptPsi2S/hMvsPt"), v12.M(), v12.Pt());
                } else {
                  fRegistry.fill(HIST("Generated/sm/PromptPsi2S/hMvsPt"), v12.M(), v12.Pt());
                }
                break;
              }
              default:
                break;
            }
          }
        } else if (hfee_type > -1) {
          auto mp1 = mcparticles.iteratorAt(t1.mothersIds()[0]);
          auto mp2 = mcparticles.iteratorAt(t2.mothersIds()[0]);
          switch (hfee_type) {
            case static_cast<int>(EM_HFeeType::kCe_Ce): {
              fRegistry.fill(HIST("Generated/ccbar/c2mu_c2mu/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Generated/ccbar/c2mu_c2mu/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Generated/ccbar/c2mu_c2mu/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/ccbar/c2mu_c2mu/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBe_Be): {
              fRegistry.fill(HIST("Generated/bbbar/b2mu_b2mu/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if (isBeautyMeson(mp1) && isBeautyMeson(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2mu_b2mu/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if (isBeautyBaryon(mp1) && isBeautyBaryon(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2mu_b2mu/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2mu_b2mu/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_BCe): {
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2c2mu/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if (isCharmMeson(mp1) && isCharmMeson(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2mu_b2mu/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if (isCharmBaryon(mp1) && isCharmBaryon(mp2)) {
                fRegistry.fill(HIST("Generated/bbbar/b2mu_b2mu/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2mu_b2mu/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): { // ULS
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_sameb/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_sameb/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_sameb/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_sameb/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
              LOGF(info, "You should not see kBCe_Be_DiffB in ULS. Good luck.");
              break;
            default:
              break;
          }
        } // end of HF evaluation
      }   // end of true ULS pair loop

      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) {
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isInAcceptance(t1) || !isInAcceptance(t2)) {
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
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.Rapidity() < minY || maxY < v12.Rapidity()) {
          continue;
        }
        if (hfee_type > -1) {
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
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            default:
              break;
          }
        }
      } // end of true LS++ pair loop

      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) {
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (!isInAcceptance(t1) || !isInAcceptance(t2)) {
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
        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.Rapidity() < minY || maxY < v12.Rapidity()) {
          continue;
        }
        if (hfee_type > -1) {
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
              fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/hadron_hadron/hMvsPt"), v12.M(), v12.Pt());
              if ((isCharmMeson(mp1) && isBeautyMeson(mp2)) || (isCharmMeson(mp2) && isBeautyMeson(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/meson_meson/hMvsPt"), v12.M(), v12.Pt());
              } else if ((isCharmBaryon(mp1) && isBeautyBaryon(mp2)) || (isCharmBaryon(mp2) && isBeautyBaryon(mp1))) {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/baryon_baryon/hMvsPt"), v12.M(), v12.Pt());
              } else {
                fRegistry.fill(HIST("Generated/bbbar/b2c2mu_b2mu_diffb/meson_baryon/hMvsPt"), v12.M(), v12.Pt());
              }
              break;
            }
            default:
              break;
          }
        }
      } // end of true LS++ pair loop

    } // end of collision loop
  }
  PROCESS_SWITCH(dimuonQCMC, processGen, "run genrated info", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(dimuonQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dimuonQCMC>(cfgc, TaskName{"dimuon-qc-mc"})};
}
