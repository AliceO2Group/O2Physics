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

#ifndef PWGEM_DILEPTON_CORE_SINGLETRACKQCMC_H_
#define PWGEM_DILEPTON_CORE_SINGLETRACKQCMC_H_

#include <map>
#include <string>
#include <vector>
#include <unordered_map>

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

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronEMEventIds, aod::EMAmbiguousElectronSelfIds, aod::EMPrimaryElectronsPrefilterBit, aod::EMPrimaryElectronMCLabels>;
using MyMCElectron = MyMCElectrons::iterator;
using FilteredMyMCElectrons = soa::Filtered<MyMCElectrons>;

using MyMCMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMAmbiguousMuonSelfIds, aod::EMPrimaryMuonMCLabels>;
using MyMCMuon = MyMCMuons::iterator;
using FilteredMyMCMuons = soa::Filtered<MyMCMuons>;

using MySmearedElectrons = soa::Join<aod::EMMCParticles, aod::SmearedElectrons>;
using MySmearedElectron = MySmearedElectrons::iterator;

using MySmearedMuons = soa::Join<aod::EMMCParticles, aod::SmearedMuons>;
using MySmearedMuon = MySmearedMuons::iterator;

template <o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType pairtype, typename TLeptons, typename TSmearedMCParticles>
struct SingleTrackQCMC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<int> cfgNtracksPV08Min{"cfgNtracksPV08Min", -1, "min. multNTracksPV"};
  Configurable<int> cfgNtracksPV08Max{"cfgNtracksPV08Max", static_cast<int>(1e+9), "max. multNTracksPV"};
  Configurable<bool> cfgFillQA{"cfgFillQA", false, "flag to fill QA histograms"};
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
    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.2, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<float> cfg_min_phi_track{"cfg_min_phi_track", 0.f, "max phi for single track"};
    Configurable<float> cfg_max_phi_track{"cfg_max_phi_track", 6.3, "max phi for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+10, "max chi2"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 1e+10, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
  } dimuoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt{"min_mcPt", 0.2, "min. MC pT"};
    Configurable<float> max_mcPt{"max_mcPt", 1e+10, "max. MC pT"};
    Configurable<float> min_mcEta{"min_mcEta", -0.8, "max. MC eta"};
    Configurable<float> max_mcEta{"max_mcEta", +0.8, "max. MC eta"};
  } mctrackcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false}; // 1 HistogramRegistry can keep up to 512 histograms
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view lepton_source_types[9] = {"lf/", "Photon/", "PromptJPsi/", "NonPromptJPsi/", "PromptPsi2S/", "NonPromptPsi2S/", "c2l/", "b2l/", "b2c2l/"};

  ~SingleTrackQCMC()
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
      const AxisSpec axis_charge_gen{3, -1.5, +1.5, "true charge"};
      std::string dca_axis_title = "DCA_{e}^{3D} (#sigma)";
      if (cfgUseDCAxy) {
        dca_axis_title = "DCA_{e}^{XY} (#sigma)";
      }
      const AxisSpec axis_dca{{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, dca_axis_title};

      // generated info
      fRegistry.add("Generated/lf/hs", "gen. single electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_charge_gen}, true);
      fRegistry.addClone("Generated/lf/", "Generated/PromptJPsi/");
      fRegistry.addClone("Generated/lf/", "Generated/NonPromptJPsi/");
      fRegistry.addClone("Generated/lf/", "Generated/PromptPsi2S/");
      fRegistry.addClone("Generated/lf/", "Generated/NonPromptPsi2S/");
      fRegistry.addClone("Generated/lf/", "Generated/c2l/");
      fRegistry.addClone("Generated/lf/", "Generated/b2l/");
      fRegistry.addClone("Generated/lf/", "Generated/b2c2l/");

      // track info
      fRegistry.add("Track/lf/positive/hs", "rec. single electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_dca, axis_charge_gen}, true);
      if (cfgFillQA) {
        fRegistry.add("Track/lf/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
        fRegistry.add("Track/lf/positive/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
        fRegistry.add("Track/lf/positive/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
        fRegistry.add("Track/lf/positive/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
        fRegistry.add("Track/lf/positive/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
        fRegistry.add("Track/lf/positive/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
        fRegistry.add("Track/lf/positive/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
        fRegistry.add("Track/lf/positive/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
        fRegistry.add("Track/lf/positive/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
        fRegistry.add("Track/lf/positive/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
        fRegistry.add("Track/lf/positive/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
        fRegistry.add("Track/lf/positive/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
        fRegistry.add("Track/lf/positive/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
        fRegistry.add("Track/lf/positive/hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{200, 0, 10}, {200, -1.0f, 1.0f}}, true);
        fRegistry.add("Track/lf/positive/hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{200, 0, 10}, {100, -0.05f, 0.05f}}, true);
        fRegistry.add("Track/lf/positive/hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{200, 0, 10}, {100, -0.05f, 0.05f}}, true);
      }
      fRegistry.addClone("Track/lf/positive/", "Track/lf/negative/");
      fRegistry.addClone("Track/lf/", "Track/Photon/"); // this is not for efficiency! only for contamination. We don't store generated photon conversions.
      fRegistry.addClone("Track/lf/", "Track/PromptJPsi/");
      fRegistry.addClone("Track/lf/", "Track/NonPromptJPsi/");
      fRegistry.addClone("Track/lf/", "Track/PromptPsi2S/");
      fRegistry.addClone("Track/lf/", "Track/NonPromptPsi2S/");
      fRegistry.addClone("Track/lf/", "Track/c2l/");
      fRegistry.addClone("Track/lf/", "Track/b2l/");
      fRegistry.addClone("Track/lf/", "Track/b2c2l/");

      if (cfgFillQA) {
        fRegistry.add("Track/PID/positive/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
        fRegistry.add("Track/PID/positive/hTOFbeta", "TOF #beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
        fRegistry.add("Track/PID/positive/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0.f, 10.f}, {160, 0, 16}}, false);
        fRegistry.add("Track/PID/positive/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTOFNsigmaMu", "TOF n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTOFNsigmaKa", "TOF n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.add("Track/PID/positive/hTOFNsigmaPr", "TOF n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
        fRegistry.addClone("Track/PID/positive/", "Track/PID/negative/");
      }
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      const AxisSpec axis_pt{ConfPtlBins, "p_{T,#mu} (GeV/c)"};
      const AxisSpec axis_eta{25, -4.5, -2.0, "#eta_{#mu}"};
      const AxisSpec axis_phi{36, 0.0, 2 * M_PI, "#varphi_{#mu} (rad.)"};
      const AxisSpec axis_dca{{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCA_{#mu}^{XY} (#sigma)"};
      const AxisSpec axis_charge_gen{3, -1.5, +1.5, "true charge"};

      // generated info
      fRegistry.add("Generated/lf/hs", "gen. single muon", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_charge_gen}, true);
      fRegistry.addClone("Generated/lf/", "Generated/PromptJPsi/");
      fRegistry.addClone("Generated/lf/", "Generated/NonPromptJPsi/");
      fRegistry.addClone("Generated/lf/", "Generated/PromptPsi2S/");
      fRegistry.addClone("Generated/lf/", "Generated/NonPromptPsi2S/");
      fRegistry.addClone("Generated/lf/", "Generated/c2l/");
      fRegistry.addClone("Generated/lf/", "Generated/b2l/");
      fRegistry.addClone("Generated/lf/", "Generated/b2c2l/");

      // track info
      fRegistry.add("Track/lf/positive/hs", "rec. single muon", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_dca, axis_charge_gen}, true);
      if (cfgFillQA) {
        fRegistry.add("Track/lf/positive/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
        fRegistry.add("Track/lf/positive/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
        fRegistry.add("Track/lf/positive/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
        fRegistry.add("Track/lf/positive/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
        fRegistry.add("Track/lf/positive/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
        fRegistry.add("Track/lf/positive/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{200, 0, 10}, {200, 0., 400}}, false);
        fRegistry.add("Track/lf/positive/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
        fRegistry.add("Track/lf/positive/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
        fRegistry.add("Track/lf/positive/hPDCA", "pDCA;p_{T} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
        fRegistry.add("Track/lf/positive/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
        fRegistry.add("Track/lf/positive/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
        fRegistry.add("Track/lf/positive/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
        fRegistry.add("Track/lf/positive/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
        fRegistry.add("Track/lf/positive/hPtGen_DeltaPtOverPtGen", "muon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{200, 0, 10}, {200, -1.0f, 1.0f}}, true);
        fRegistry.add("Track/lf/positive/hPtGen_DeltaEta", "muon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{200, 0, 10}, {100, -0.05f, 0.05f}}, true);
        fRegistry.add("Track/lf/positive/hPtGen_DeltaPhi", "muon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{200, 0, 10}, {100, -0.05f, 0.05f}}, true);
      }
      fRegistry.addClone("Track/lf/positive/", "Track/lf/negative/");
      fRegistry.addClone("Track/lf/", "Track/Photon/"); // this is not for efficiency! only for contamination. We don't store generated photon conversions.
      fRegistry.addClone("Track/lf/", "Track/PromptJPsi/");
      fRegistry.addClone("Track/lf/", "Track/NonPromptJPsi/");
      fRegistry.addClone("Track/lf/", "Track/PromptPsi2S/");
      fRegistry.addClone("Track/lf/", "Track/NonPromptPsi2S/");
      fRegistry.addClone("Track/lf/", "Track/c2l/");
      fRegistry.addClone("Track/lf/", "Track/b2l/");
      fRegistry.addClone("Track/lf/", "Track/b2c2l/");
    }
  }

  int pdg_lepton = 0;
  void init(InitContext&)
  {
    if (doprocessQCMC && doprocessQCMC_Smeared) {
      LOGF(fatal, "Cannot enable processQCMC and processQCMC_Smeared at the same time. Please choose one.");
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    DefineEMEventCut();
    addhistograms();

    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      pdg_lepton = 11;
      DefineDielectronCut();
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      pdg_lepton = 13;
      DefineDimuonCut();
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

    if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
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

  template <bool isSmeared, typename T>
  bool isInAcceptance(T const& lepton)
  {
    float pt = 0.f, eta = 0.f;
    if constexpr (isSmeared) {
      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        pt = lepton.ptSmeared();
        eta = lepton.etaSmeared();
      } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
        if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          pt = lepton.ptSmeared_sa_muon();
          eta = lepton.etaSmeared_sa_muon();
        } else if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          pt = lepton.ptSmeared_gl_muon();
          eta = lepton.etaSmeared_gl_muon();
        } else {
          pt = lepton.pt();
          eta = lepton.eta();
        }
      }
    } else {
      pt = lepton.pt();
      eta = lepton.eta();
    }

    if ((mctrackcuts.min_mcPt < pt && pt < mctrackcuts.max_mcPt) && (mctrackcuts.min_mcEta < eta && eta < mctrackcuts.max_mcEta)) {
      return true;
    } else {
      return false;
    }
  }

  template <int lepton_source_id, typename TMCParticles, typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      fillElectronInfo<lepton_source_id, TMCParticles>(track);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      fillMuonInfo<lepton_source_id, TMCParticles>(track);
    }
  }

  template <int lepton_source_id, typename TMCParticles, typename TTrack>
  void fillElectronInfo(TTrack const& track)
  {
    auto mctrack = track.template emmcparticle_as<TMCParticles>();
    float dca = dca3DinSigma(track);
    if (cfgUseDCAxy) {
      dca = abs(track.dcaXY() / std::sqrt(track.cYY()));
    }

    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[track.globalIndex()];
    }
    // LOGF(info, "map_weight[%d] = %f", track.globalIndex(), weight);

    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hs"), track.pt(), track.eta(), track.phi(), dca, -mctrack.pdgCode() / pdg_lepton, weight);
      if (cfgFillQA) {
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hQoverPt"), track.sign() / track.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAxyz"), track.dcaXY(), track.dcaZ());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hNclsITS"), track.itsNCls());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hNclsTPC"), track.tpcNClsFound());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hNcrTPC"), track.tpcNClsCrossedRows());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hChi2TPC"), track.tpcChi2NCl());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hChi2ITS"), track.itsChi2NCl());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hITSClusterMap"), track.itsClusterMap());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());

        fRegistry.fill(HIST("Track/PID/positive/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
        fRegistry.fill(HIST("Track/PID/positive/hTOFbeta"), track.p(), track.beta());
        fRegistry.fill(HIST("Track/PID/positive/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
        fRegistry.fill(HIST("Track/PID/positive/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        fRegistry.fill(HIST("Track/PID/positive/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
        fRegistry.fill(HIST("Track/PID/positive/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        fRegistry.fill(HIST("Track/PID/positive/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
        fRegistry.fill(HIST("Track/PID/positive/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        fRegistry.fill(HIST("Track/PID/positive/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
        fRegistry.fill(HIST("Track/PID/positive/hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
        fRegistry.fill(HIST("Track/PID/positive/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
        fRegistry.fill(HIST("Track/PID/positive/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
        fRegistry.fill(HIST("Track/PID/positive/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
      }
    } else {
      fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hs"), track.pt(), track.eta(), track.phi(), dca, -mctrack.pdgCode() / pdg_lepton, weight);
      if (cfgFillQA) {
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hQoverPt"), track.sign() / track.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAxyz"), track.dcaXY(), track.dcaZ());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAxyRes_Pt"), track.pt(), sqrt(track.cYY()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAzRes_Pt"), track.pt(), sqrt(track.cZZ()) * 1e+4);  // convert cm to um
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hNclsITS"), track.itsNCls());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hNclsTPC"), track.tpcNClsFound());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hNcrTPC"), track.tpcNClsCrossedRows());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hChi2TPC"), track.tpcChi2NCl());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hChi2ITS"), track.itsChi2NCl());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hITSClusterMap"), track.itsClusterMap());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());

        fRegistry.fill(HIST("Track/PID/negative/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
        fRegistry.fill(HIST("Track/PID/negative/hTOFbeta"), track.p(), track.beta());
        fRegistry.fill(HIST("Track/PID/negative/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
        fRegistry.fill(HIST("Track/PID/negative/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        fRegistry.fill(HIST("Track/PID/negative/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
        fRegistry.fill(HIST("Track/PID/negative/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        fRegistry.fill(HIST("Track/PID/negative/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
        fRegistry.fill(HIST("Track/PID/negative/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        fRegistry.fill(HIST("Track/PID/negative/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
        fRegistry.fill(HIST("Track/PID/negative/hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
        fRegistry.fill(HIST("Track/PID/negative/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
        fRegistry.fill(HIST("Track/PID/negative/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
        fRegistry.fill(HIST("Track/PID/negative/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
      }
    }
  }

  template <int lepton_source_id, typename TMCParticles, typename TTrack>
  void fillMuonInfo(TTrack const& track)
  {
    auto mctrack = track.template emmcparticle_as<TMCParticles>();
    float dca_xy = fwdDcaXYinSigma(track);

    float weight = 1.f;
    if (cfgApplyWeightTTCA) {
      weight = map_weight[track.globalIndex()];
    }
    // LOGF(info, "map_weight[%d] = %f", track.globalIndex(), weight);

    if (track.sign() > 0) {
      fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hs"), track.pt(), track.eta(), track.phi(), dca_xy, -mctrack.pdgCode() / pdg_lepton, weight);
      if (cfgFillQA) {
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hQoverPt"), track.sign() / track.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hTrackType"), track.trackType());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXX()), track.fwdDcaY() / std::sqrt(track.cYY()));
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXX()) * 1e+4);
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4);
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hNclsMCH"), track.nClusters());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hNclsMFT"), track.nClustersMFT());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hPDCA"), track.pt(), track.pDca());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hChi2"), track.chi2());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hMFTClusterMap"), track.mftClusterMap());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("positive/hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());
      }
    } else {
      fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hs"), track.pt(), track.eta(), track.phi(), dca_xy, -mctrack.pdgCode() / pdg_lepton, weight);
      if (cfgFillQA) {
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hQoverPt"), track.sign() / track.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hTrackType"), track.trackType());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXX()), track.fwdDcaY() / std::sqrt(track.cYY()));
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXX()) * 1e+4);
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4);
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hNclsMCH"), track.nClusters());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hNclsMFT"), track.nClustersMFT());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hPDCA"), track.pt(), track.pDca());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hChi2"), track.chi2());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hChi2MatchMCHMID"), track.chi2MatchMCHMID());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hMFTClusterMap"), track.mftClusterMap());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
        fRegistry.fill(HIST("Track/") + HIST(lepton_source_types[lepton_source_id]) + HIST("negative/hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());
      }
    }
  }

  template <typename TCollisions, typename TTracks, typename TPreslice, typename TCut, typename TMCCollisions, typename TMCParticles>
  void runQCMC(TCollisions const& collisions, TTracks const& tracks, TPreslice const& perCollision, TCut const& cut, TMCCollisions const&, TMCParticles const& mcparticles)
  {
    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      // auto mccollision = collision.template emmcevent_as<TMCCollisions>();
      // if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
      //   continue;
      // }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());

      for (auto& track : tracks_per_coll) {
        auto mctrack = track.template emmcparticle_as<TMCParticles>();
        if (abs(mctrack.pdgCode()) != pdg_lepton) {
          continue;
        }

        auto mccollision_from_track = mctrack.template emmcevent_as<TMCCollisions>();
        if (cfgEventGeneratorType >= 0 && mccollision_from_track.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }

        // if (mctrack.emmceventId() != collision.emmceventId()) {
        //   continue;
        // }

        if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
          if (dielectroncuts.cfg_pid_scheme == static_cast<int>(DielectronCut::PIDSchemes::kPIDML)) {
            if (!cut.template IsSelectedTrack<true>(track, collision)) {
              continue;
            }
          } else { // cut-based
            if (!cut.template IsSelectedTrack(track)) {
              continue;
            }
          }
        } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
          if (!cut.template IsSelectedTrack(track)) {
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
            fillTrackInfo<0, TMCParticles>(track);
          } else if (pdg_mother == 443) {
            if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
              fillTrackInfo<3, TMCParticles>(track);
            } else {
              fillTrackInfo<2, TMCParticles>(track);
            }
          } else if (pdg_mother == 100443) {
            if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
              fillTrackInfo<5, TMCParticles>(track);
            } else {
              fillTrackInfo<4, TMCParticles>(track);
            }
          } else if (IsFromBeauty(mctrack, mcparticles) > 0) { // b is found in full decay chain.
            if (IsFromCharm(mctrack, mcparticles) > 0) {       // c is found in full decay chain.
              fillTrackInfo<8, TMCParticles>(track);
            } else {
              fillTrackInfo<7, TMCParticles>(track);
            }
          } else if (IsFromCharm(mctrack, mcparticles) > 0) { // c is found in full decay chain. Not from b.
            fillTrackInfo<6, TMCParticles>(track);
          }
        } else {
          fillTrackInfo<1, TMCParticles>(track);
        }
      } // end of track loop

    } // end of collision loop

  } // end of process

  template <bool isSmeared, typename TCollisions, typename TMCLeptons, typename TMCCollisions, typename TMCParticles>
  void runGenInfo(TCollisions const& collisions, TMCLeptons const& leptonsMC, TMCCollisions const&, TMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      float centralities[4] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      auto mccollision = collision.template emmcevent_as<TMCCollisions>();
      // LOGF(info, "mccollision.getGeneratorId() = %d", mccollision.getGeneratorId());
      // LOGF(info, "mccollision.getSubGeneratorId() = %d", mccollision.getSubGeneratorId());
      // LOGF(info, "mccollision.getSourceId() = %d", mccollision.getSourceId());
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted

      auto leptonsMC_per_coll = leptonsMC.sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);

      for (auto& lepton : leptonsMC_per_coll) {
        if (!(lepton.isPhysicalPrimary() || lepton.producedByGenerator())) {
          continue;
        }
        if (!isInAcceptance<isSmeared>(lepton)) {
          continue;
        }
        if (!lepton.has_mothers()) {
          continue;
        }
        auto mcmother = mcparticles.iteratorAt(lepton.mothersIds()[0]);
        int pdg_mother = abs(mcmother.pdgCode());

        float pt = 0.f, eta = 0.f, phi = 0.f;
        if constexpr (isSmeared) {
          if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
            pt = lepton.ptSmeared();
            eta = lepton.etaSmeared();
            phi = lepton.phiSmeared();
          } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
            if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
              pt = lepton.ptSmeared_sa_muon();
              eta = lepton.etaSmeared_sa_muon();
              phi = lepton.phiSmeared_sa_muon();
            } else if (dimuoncuts.cfg_track_type == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
              pt = lepton.ptSmeared_gl_muon();
              eta = lepton.etaSmeared_gl_muon();
              phi = lepton.phiSmeared_gl_muon();
            } else {
              pt = lepton.pt();
              eta = lepton.eta();
              phi = lepton.phi();
            }
          }
        } else {
          pt = lepton.pt();
          eta = lepton.eta();
          phi = lepton.phi();
        }

        if (pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 331 || pdg_mother == 113 || pdg_mother == 223 || pdg_mother == 333) {
          fRegistry.fill(HIST("Generated/lf/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
        } else if (pdg_mother == 443) {
          if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
            fRegistry.fill(HIST("Generated/NonPromptJPsi/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
          } else {
            fRegistry.fill(HIST("Generated/PromptJPsi/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
          }
        } else if (pdg_mother == 100443) {
          if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
            fRegistry.fill(HIST("Generated/NonPromptPsi2S/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
          } else {
            fRegistry.fill(HIST("Generated/PromptPsi2S/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
          }
        } else if (IsFromBeauty(lepton, mcparticles) > 0) { // b is found in full decay chain.
          if (IsFromCharm(lepton, mcparticles) > 0) {       // c is found in full decay chain.
            fRegistry.fill(HIST("Generated/b2c2l/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
          } else {
            fRegistry.fill(HIST("Generated/b2l/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
          }
        } else if (IsFromCharm(lepton, mcparticles) > 0) { // c is found in full decay chain. Not from b.
          fRegistry.fill(HIST("Generated/c2l/hs"), pt, eta, phi, -lepton.pdgCode() / pdg_lepton);
        }
      }

    } // end of collision loop
  }

  std::unordered_map<int, float> map_weight; // map of track global index -> weight
  template <typename TCollisions, typename TTracks, typename TPreslice, typename TCut, typename TMCCollisions, typename TMCParticles>
  void fillTrackWeightMap(TCollisions const& collisions, TTracks const& tracks, TPreslice const& perCollision, TCut const& cut, TMCCollisions const&, TMCParticles const&)
  {
    std::vector<int> passed_trackIds;
    passed_trackIds.reserve(tracks.size());
    for (auto& collision : collisions) {
      float centralities[4] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      // auto mccollision = collision.template emmcevent_as<TMCCollisions>();
      // if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
      //   continue;
      // }

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());

      if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
        for (auto& track : tracks_per_coll) {
          auto mctrack = track.template emmcparticle_as<TMCParticles>();
          auto mccollision_from_track = mctrack.template emmcevent_as<TMCCollisions>();
          if (cfgEventGeneratorType >= 0 && mccollision_from_track.getSubGeneratorId() != cfgEventGeneratorType) {
            continue;
          }

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
          auto mctrack = track.template emmcparticle_as<TMCParticles>();
          auto mccollision_from_track = mctrack.template emmcevent_as<TMCCollisions>();
          if (cfgEventGeneratorType >= 0 && mccollision_from_track.getSubGeneratorId() != cfgEventGeneratorType) {
            continue;
          }

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
  Preslice<MyMCElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Filter trackFilter_electron = dielectroncuts.cfg_min_pt_track < o2::aod::track::pt && dielectroncuts.cfg_min_eta_track < o2::aod::track::eta && o2::aod::track::eta < dielectroncuts.cfg_max_eta_track && dielectroncuts.cfg_min_phi_track < o2::aod::track::phi && o2::aod::track::phi < dielectroncuts.cfg_max_phi_track && o2::aod::track::tpcChi2NCl < dielectroncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dielectroncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dielectroncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dielectroncuts.cfg_max_dcaz;
  Filter pidFilter_electron = (dielectroncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dielectroncuts.cfg_max_TPCNsigmaEl) && (o2::aod::pidtpc::tpcNSigmaPi < dielectroncuts.cfg_min_TPCNsigmaPi || dielectroncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi);
  Filter ttcaFilter_electron = ifnode(dielectroncuts.enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);

  Preslice<MyMCMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Filter trackFilter_muon = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track && dimuoncuts.cfg_min_phi_track < o2::aod::fwdtrack::phi && o2::aod::fwdtrack::phi < dimuoncuts.cfg_max_phi_track;
  Filter ttcaFilter_muon = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_multiplicity = cfgNtracksPV08Min <= o2::aod::mult::multNTracksPV && o2::aod::mult::multNTracksPV < cfgNtracksPV08Max;
  Filter collisionFilter_occupancy = eventcuts.cfgOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgOccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Partition<aod::EMMCParticles> electronsMC = nabs(o2::aod::mcparticle::pdgCode) == 11; // e+, e-
  Partition<aod::EMMCParticles> muonsMC = nabs(o2::aod::mcparticle::pdgCode) == 13;     // mu+, mu-
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;

  void processQCMC(FilteredMyCollisions const& collisions, aod::EMMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles, TLeptons const& tracks)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap(collisions, tracks, perCollision_electron, fDielectronCut, mccollisions, mcparticles);
      }
      runQCMC(collisions, tracks, perCollision_electron, fDielectronCut, mccollisions, mcparticles);
      runGenInfo<false>(collisions, electronsMC, mccollisions, mcparticles);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap(collisions, tracks, perCollision_muon, fDimuonCut, mccollisions, mcparticles);
      }
      runQCMC(collisions, tracks, perCollision_muon, fDimuonCut, mccollisions, mcparticles);
      runGenInfo<false>(collisions, muonsMC, mccollisions, mcparticles);
    }
    map_weight.clear();
  }
  PROCESS_SWITCH(SingleTrackQCMC, processQCMC, "run single track QC MC", true);

  Partition<MySmearedElectrons> electronsMC_smeared = nabs(o2::aod::mcparticle::pdgCode) == 11; // e+, e-
  Partition<MySmearedMuons> muonsMC_smeared = nabs(o2::aod::mcparticle::pdgCode) == 13;         // mu+, mu-

  void processQCMC_Smeared(FilteredMyCollisions const& collisions, aod::EMMCEvents const& mccollisions, TLeptons const& tracks, TSmearedMCParticles const& mcparticles_smeared)
  {
    if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDielectron) {
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap(collisions, tracks, perCollision_electron, fDielectronCut, mccollisions, mcparticles_smeared);
      }
      runQCMC(collisions, tracks, perCollision_electron, fDielectronCut, mccollisions, mcparticles_smeared);
      runGenInfo<true>(collisions, electronsMC_smeared, mccollisions, mcparticles_smeared);
    } else if constexpr (pairtype == o2::aod::pwgem::dilepton::utils::pairutil::DileptonPairType::kDimuon) {
      if (cfgApplyWeightTTCA) {
        fillTrackWeightMap(collisions, tracks, perCollision_muon, fDimuonCut, mccollisions, mcparticles_smeared);
      }
      runQCMC(collisions, tracks, perCollision_muon, fDimuonCut, mccollisions, mcparticles_smeared);
      runGenInfo<true>(collisions, muonsMC_smeared, mccollisions, mcparticles_smeared);
    }
    map_weight.clear();
  }
  PROCESS_SWITCH(SingleTrackQCMC, processQCMC_Smeared, "run single track QC MC with smearing", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(SingleTrackQCMC, processDummy, "Dummy function", false);
};

#endif // PWGEM_DILEPTON_CORE_SINGLETRACKQCMC_H_
