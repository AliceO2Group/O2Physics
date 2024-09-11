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
// This code runs loop over electrons for QC.
//    Please write to: daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_CORE_SINGLETRACKQC_H_
#define PWGEM_DILEPTON_CORE_SINGLETRACKQC_H_

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

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
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsWithSWT = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMSWTriggerInfos>;
using MyCollisionWithSWT = MyCollisionsWithSWT::iterator;

using MyElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyElectron = MyElectrons::iterator;
using FilteredMyElectrons = soa::Filtered<MyElectrons>;
using FilteredMyElectron = FilteredMyElectrons::iterator;

using MyMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds>;
using MyMuon = MyMuons::iterator;
using FilteredMyMuons = soa::Filtered<MyMuons>;
using FilteredMyMuon = FilteredMyMuons::iterator;

template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename... Types>
struct SingleTrackQC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2, NTPV:3"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<std::string> cfg_swt_name{"cfg_swt_name", "fHighTrackMult", "desired software trigger name"}; // 1 trigger name per 1 task. fHighTrackMult, fHighFt0Mult
  Configurable<int> cfgNtracksPV08Min{"cfgNtracksPV08Min", -1, "min. multNTracksPV"};
  Configurable<int> cfgNtracksPV08Max{"cfgNtracksPV08Max", static_cast<int>(1e+9), "max. multNTracksPV"};
  Configurable<bool> cfgApplyWeightTTCA{"cfgApplyWeightTTCA", false, "flag to apply weighting by 1/N"};
  Configurable<bool> cfgUseDCAxy{"cfgUseDCAxy", false, "flag to use DCAxy, instead of DCA3D"};

  ConfigurableAxis ConfPtlBins{"ConfPtlBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTl bins for output histograms"};

  EMEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -1, "min. occupancy"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
  } eventcuts;

  DielectronCut fDielectronCut;
  struct : ConfigurableGroup {
    std::string prefix = "dielectroncut_group";

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "max eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "max phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 100, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};

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
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "max phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
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

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false}; // 1 HistogramRegistry can keep up to 512 histograms
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};

  ~SingleTrackQC()
  {
    if (eid_bdt) {
      delete eid_bdt;
    }
  }

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms<-1>(&fRegistry);

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      const AxisSpec axis_pt{ConfPtlBins, "p_{T,e} (GeV/c)"};
      const AxisSpec axis_eta{20, -1.0, +1.0, "#eta_{e}"};
      const AxisSpec axis_phi{36, 0.0, 2 * M_PI, "#varphi_{e} (rad.)"};
      std::string dca_axis_title = "DCA_{e}^{3D} (#sigma)";
      if (cfgUseDCAxy) {
        dca_axis_title = "DCA_{e}^{XY} (#sigma)";
      }
      const AxisSpec axis_dca{{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, dca_axis_title};

      // track info
      fRegistry.add("Track/positive/hs", "rec. single electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_dca}, true);
      fRegistry.add("Track/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
      fRegistry.add("Track/positive/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/positive/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/positive/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
      fRegistry.add("Track/positive/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
      fRegistry.add("Track/positive/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/positive/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/positive/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/positive/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/positive/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/positive/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("Track/positive/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/positive/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("Track/positive/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/positive/hTOFbeta", "TOF #beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
      fRegistry.add("Track/positive/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda);", kTH2F, {{1000, 0.f, 10.f}, {160, 0, 16}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTOFNsigmaMu", "TOF n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTOFNsigmaKa", "TOF n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/positive/hTOFNsigmaPr", "TOF n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.addClone("Track/positive/", "Track/negative/");
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      const AxisSpec axis_pt{ConfPtlBins, "p_{T,#mu} (GeV/c)"};
      const AxisSpec axis_eta{25, -4.5, -2.0, "#eta_{#mu}"};
      const AxisSpec axis_phi{36, 0.0, 2 * M_PI, "#varphi_{#mu} (rad.)"};
      const AxisSpec axis_dca{{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCA_{#mu}^{XY} (#sigma)"};

      // track info
      fRegistry.add("Track/positive/hs", "rec. single muon", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_dca}, true);
      fRegistry.add("Track/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
      fRegistry.add("Track/positive/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
      fRegistry.add("Track/positive/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/positive/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/positive/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
      fRegistry.add("Track/positive/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
      fRegistry.add("Track/positive/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
      fRegistry.add("Track/positive/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
      fRegistry.add("Track/positive/hPDCA", "pDCA;p_{T} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
      fRegistry.add("Track/positive/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/positive/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/positive/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/positive/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
      fRegistry.addClone("Track/positive/", "Track/negative/");
    }
  }

  int mRunNumber;
  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    DefineEMEventCut();
    DefineDielectronCut();
    DefineDimuonCut();
    addhistograms();
    mRunNumber = 0;

    if (doprocessQC_TriggeredData) {
      fRegistry.add("Event/hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
    }
  }

  template <bool isTriggerAnalysis, typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    mRunNumber = collision.runNumber();

    if constexpr (isTriggerAnalysis) {
      LOGF(info, "Trigger analysis is enabled. Desired trigger name = %s", cfg_swt_name.value);
      LOGF(info, "total inspected TVX events = %d in run number %d", collision.nInspectedTVX(), collision.runNumber());
      fRegistry.fill(HIST("Event/hNInspectedTVX"), collision.runNumber(), collision.nInspectedTVX());
    }
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
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
  }

  o2::ml::OnnxModel* eid_bdt = nullptr;
  void DefineDielectronCut()
  {
    fDielectronCut = DielectronCut("fDielectronCut", "fDielectronCut");

    // for track
    fDielectronCut.SetTrackPtRange(dielectroncuts.cfg_min_pt_track, 1e+10f);
    fDielectronCut.SetTrackEtaRange(dielectroncuts.cfg_min_eta_track, dielectroncuts.cfg_max_eta_track);
    fDielectronCut.SetTrackPhiRange(dielectroncuts.cfg_min_phi_track, dielectroncuts.cfg_max_phi_track);
    fDielectronCut.SetMinNClustersTPC(dielectroncuts.cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(dielectroncuts.cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetChi2PerClusterTPC(0.0, dielectroncuts.cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, dielectroncuts.cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(dielectroncuts.cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITS(0, 16);
    fDielectronCut.SetMaxDcaXY(dielectroncuts.cfg_max_dcaxy);
    fDielectronCut.SetMaxDcaZ(dielectroncuts.cfg_max_dcaz);
    fDielectronCut.RequireITSibAny(dielectroncuts.cfg_require_itsib_any);
    fDielectronCut.RequireITSib1st(dielectroncuts.cfg_require_itsib_1st);

    // for eID
    fDielectronCut.SetPIDScheme(dielectroncuts.cfg_pid_scheme);
    fDielectronCut.SetTPCNsigmaElRange(dielectroncuts.cfg_min_TPCNsigmaEl, dielectroncuts.cfg_max_TPCNsigmaEl);
    fDielectronCut.SetTPCNsigmaMuRange(dielectroncuts.cfg_min_TPCNsigmaMu, dielectroncuts.cfg_max_TPCNsigmaMu);
    fDielectronCut.SetTPCNsigmaPiRange(dielectroncuts.cfg_min_TPCNsigmaPi, dielectroncuts.cfg_max_TPCNsigmaPi);
    fDielectronCut.SetTPCNsigmaKaRange(dielectroncuts.cfg_min_TPCNsigmaKa, dielectroncuts.cfg_max_TPCNsigmaKa);
    fDielectronCut.SetTPCNsigmaPrRange(dielectroncuts.cfg_min_TPCNsigmaPr, dielectroncuts.cfg_max_TPCNsigmaPr);
    fDielectronCut.SetTOFNsigmaElRange(dielectroncuts.cfg_min_TOFNsigmaEl, dielectroncuts.cfg_max_TOFNsigmaEl);

    if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDielectronCut
      eid_bdt = new o2::ml::OnnxModel();
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
    fDimuonCut.SetTrackPhiRange(dimuoncuts.cfg_min_phi_track, dimuoncuts.cfg_max_phi_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 16);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
  }

  template <typename TTrack>
  void fillElectronInfo(TTrack const& track)
  {
    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[track.globalIndex()];
    }

    float dca = dca3DinSigma(track);
    if (cfgUseDCAxy) {
      dca = abs(track.dcaXY() / std::sqrt(track.cYY()));
    }

    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/positive/hs"), track.pt(), track.eta(), track.phi(), dca, weight);
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
      fRegistry.fill(HIST("Track/positive/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/positive/hITSClusterMap"), track.itsClusterMap());

      fRegistry.fill(HIST("Track/positive/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/positive/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/positive/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/positive/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/positive/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
    } else {
      fRegistry.fill(HIST("Track/negative/hs"), track.pt(), track.eta(), track.phi(), dca, weight);
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
      fRegistry.fill(HIST("Track/negative/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/negative/hITSClusterMap"), track.itsClusterMap());

      fRegistry.fill(HIST("Track/negative/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/negative/hTOFbeta"), track.p(), track.beta());
      fRegistry.fill(HIST("Track/negative/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
      fRegistry.fill(HIST("Track/negative/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
      fRegistry.fill(HIST("Track/negative/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
    }
  }

  template <typename TTrack>
  void fillMuonInfo(TTrack const& track)
  {
    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[track.globalIndex()];
    }
    float dca_xy = fwdDcaXYinSigma(track);
    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/positive/hs"), track.pt(), track.eta(), track.phi(), dca_xy, weight);
      fRegistry.fill(HIST("Track/positive/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/positive/hTrackType"), track.trackType());
      fRegistry.fill(HIST("Track/positive/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
      fRegistry.fill(HIST("Track/positive/hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXX()), track.fwdDcaY() / std::sqrt(track.cYY()));
      fRegistry.fill(HIST("Track/positive/hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXX()) * 1e+4);
      fRegistry.fill(HIST("Track/positive/hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4);
      fRegistry.fill(HIST("Track/positive/hNclsMCH"), track.nClusters());
      fRegistry.fill(HIST("Track/positive/hNclsMFT"), track.nClustersMFT());
      fRegistry.fill(HIST("Track/positive/hPDCA"), track.pt(), track.pDca());
      fRegistry.fill(HIST("Track/positive/hChi2"), track.chi2());
      fRegistry.fill(HIST("Track/positive/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
      fRegistry.fill(HIST("Track/positive/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
      fRegistry.fill(HIST("Track/positive/hMFTClusterMap"), track.mftClusterMap());
    } else {
      fRegistry.fill(HIST("Track/negative/hs"), track.pt(), track.eta(), track.phi(), dca_xy, weight);
      fRegistry.fill(HIST("Track/negative/hQoverPt"), track.sign() / track.pt());
      fRegistry.fill(HIST("Track/negative/hTrackType"), track.trackType());
      fRegistry.fill(HIST("Track/negative/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
      fRegistry.fill(HIST("Track/negative/hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXX()), track.fwdDcaY() / std::sqrt(track.cYY()));
      fRegistry.fill(HIST("Track/negative/hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXX()) * 1e+4);
      fRegistry.fill(HIST("Track/negative/hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4);
      fRegistry.fill(HIST("Track/negative/hNclsMCH"), track.nClusters());
      fRegistry.fill(HIST("Track/negative/hNclsMFT"), track.nClustersMFT());
      fRegistry.fill(HIST("Track/negative/hPDCA"), track.pt(), track.pDca());
      fRegistry.fill(HIST("Track/negative/hChi2"), track.chi2());
      fRegistry.fill(HIST("Track/negative/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
      fRegistry.fill(HIST("Track/negative/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
      fRegistry.fill(HIST("Track/negative/hMFTClusterMap"), track.mftClusterMap());
    }
  }

  template <bool isTriggerAnalysis, typename TCollisions, typename TTracks, typename TPreslice, typename TCut>
  void runQC(TCollisions const& collisions, TTracks const& tracks, TPreslice const& perCollision, TCut const& cut)
  {
    for (auto& collision : collisions) {
      initCCDB<isTriggerAnalysis>(collision);
      float centralities[4] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      if constexpr (isTriggerAnalysis) {
        if (!collision.swtalias_bit(o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value))) {
          continue;
        }
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        for (auto& track : tracks_per_coll) {
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut.template IsSelectedTrack<true>(track, collision)) {
              continue;
            }
          } else { // cut-based
            if (!cut.template IsSelectedTrack(track)) {
              continue;
            }
          }
          fillElectronInfo(track);
        } // end of track loop
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        for (auto& track : tracks_per_coll) {
          if (!cut.template IsSelectedTrack(track)) {
            continue;
          }
          fillMuonInfo(track);
        } // end of track loop
      }
    } // end of collision loop
  }

  std::unordered_map<int, float> map_weight;
  template <bool isTriggerAnalysis, typename TCollisions, typename TTracks, typename TPreslice, typename TCut>
  void fillTrackWeightMap(TCollisions const& collisions, TTracks const& tracks, TPreslice const& perCollision, TCut const& cut)
  {
    std::vector<int> passed_trackIds;
    passed_trackIds.reserve(tracks.size());
    for (auto& collision : collisions) {
      initCCDB<isTriggerAnalysis>(collision);
      float centralities[4] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      if constexpr (isTriggerAnalysis) {
        if (!collision.swtalias_bit(o2::aod::pwgem::dilepton::swt::aliasLabels.at(cfg_swt_name.value))) {
          continue;
        }
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        for (auto& track : tracks_per_coll) {
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut.template IsSelectedTrack<true>(track, collision)) {
              continue;
            }
          } else { // cut-based
            if (!cut.template IsSelectedTrack(track)) {
              continue;
            }
          }
          passed_trackIds.emplace_back(track.globalIndex());
        } // end of track loop
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        for (auto& track : tracks_per_coll) {
          if (!cut.template IsSelectedTrack(track)) {
            continue;
          }
          passed_trackIds.emplace_back(track.globalIndex());
        } // end of track loop
      }
    } // end of collision loop

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      for (auto& trackId : passed_trackIds) {
        auto track = tracks.rawIteratorAt(trackId);
        auto ambIds = track.ambiguousElectronsIds();
        float n = 1.f; // include myself.
        for (auto& ambId : ambIds) {
          if (std::find(passed_trackIds.begin(), passed_trackIds.end(), ambId) != passed_trackIds.end()) {
            n += 1.f;
          }
        }
        map_weight[trackId] = 1.f / n;
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      for (auto& trackId : passed_trackIds) {
        auto track = tracks.rawIteratorAt(trackId);
        auto ambIds = track.ambiguousMuonsIds();
        float n = 1.f; // include myself.
        for (auto& ambId : ambIds) {
          if (std::find(passed_trackIds.begin(), passed_trackIds.end(), ambId) != passed_trackIds.end()) {
            n += 1.f;
          }
        }
        map_weight[trackId] = 1.f / n;
      }
    }

    passed_trackIds.clear();
    passed_trackIds.shrink_to_fit();
  }

  SliceCache cache;
  Preslice<MyElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Filter trackFilter_electron = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && dielectroncuts.cfg_min_phi_track < o2::aod::track::phi && o2::aod::track::phi < dielectroncuts.cfg_max_phi_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter_electron = (dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < dielectroncuts.cfg_min_TPCNsigmaPi || dielectroncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi);
  Filter ttcaFilter_electron = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);

  Preslice<MyMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track && dimuoncuts.cfg_min_phi_track < o2::aod::fwdtrack::phi && o2::aod::fwdtrack::phi < dimuoncuts.cfg_max_phi_track;
  Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_multiplicity = cfgNtracksPV08Min <= o2::aod::mult::multNTracksPV && o2::aod::mult::multNTracksPV < cfgNtracksPV08Max;
  Filter collisionFilter_occupancy = eventcuts.cfgOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgOccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processQC(FilteredMyCollisions const& collisions, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<false>(collisions, electrons, perCollision_electron, fDielectronCut);
      }
      runQC<false>(collisions, electrons, perCollision_electron, fDielectronCut);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<false>(collisions, muons, perCollision_muon, fDimuonCut);
      }
      runQC<false>(collisions, muons, perCollision_muon, fDimuonCut);
    }

    map_weight.clear();
  }
  PROCESS_SWITCH(SingleTrackQC, processQC, "run single track QC", true);

  using FilteredMyCollisionsWithSWT = soa::Filtered<MyCollisionsWithSWT>;
  void processQC_TriggeredData(FilteredMyCollisionsWithSWT const& collisions, Types const&... args)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      auto electrons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<true>(collisions, electrons, perCollision_electron, fDielectronCut);
      }
      runQC<true>(collisions, electrons, perCollision_electron, fDielectronCut);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      auto muons = std::get<0>(std::tie(args...));
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap<true>(collisions, muons, perCollision_muon, fDimuonCut);
      }
      runQC<true>(collisions, muons, perCollision_muon, fDimuonCut);
    }

    map_weight.clear();
  }
  PROCESS_SWITCH(SingleTrackQC, processQC_TriggeredData, "run single track QC on triggered data", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(SingleTrackQC, processDummy, "Dummy function", false);
};
#endif // PWGEM_DILEPTON_CORE_SINGLETRACKQC_H_
