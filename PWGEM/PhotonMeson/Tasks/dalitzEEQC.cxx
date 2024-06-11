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
#include "Common/Core/RecoDecay.h"

#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMEventCut.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;

struct DalitzEEQC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};

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

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.5, "max mass"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.9, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3]"};
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

    // CCDB configuration for PID ML
    Configurable<std::string> BDTLocalPathGamma{"BDTLocalPathGamma", "pid_ml_xgboost.onnx", "Path to the local .onnx file"};
    Configurable<std::string> BDTPathCCDB{"BDTPathCCDB", "Users/d/dsekihat/pwgem/pidml/", "Path on CCDB"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dileptoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  bool cfgDoFlow = false;

  void init(InitContext& /*context*/)
  {
    DefineEMEventCut();
    DefineDileptonCut();
    addhistograms();
  }

  ~DalitzEEQC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry, cfgDoFlow);

    // pair info
    std::vector<double> ptbins;
    std::vector<double> massbins;

    for (int i = 0; i < 51; i++) {
      massbins.emplace_back(0.01 * (i - 0) + 0.0); // every 0.01 GeV/c2 from 0.0 to 0.5 GeV/c2
    }
    const AxisSpec axis_mass{massbins, "m_{ee} (GeV/c^{2})"};

    for (int i = 0; i < 50; i++) {
      ptbins.emplace_back(0.1 * (i - 0) + 0.0); // every 0.1 GeV/c from 0.0 to 5.0 GeV/c
    }
    for (int i = 50; i < 61; i++) {
      ptbins.emplace_back(0.5 * (i - 50) + 5.0); // every 0.5 GeV/c from 5.0 to 10 GeV/c
    }
    const AxisSpec axis_pt{ptbins, "p_{T,ee} (GeV/c)"};

    fRegistry.add("Pair/same/uls/hMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry.add("Pair/same/uls/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lspp/");
    fRegistry.addClone("Pair/same/uls/", "Pair/same/lsmm/");

    fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {100, 0., 1000}}, false);
    fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {100, 0., 1000}}, false);
    fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTOFbeta", "TOF beta;p_{in} (GeV/c);TOF #beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
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
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mass, dileptoncuts.cfg_max_mass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.SetPairDCARange(dileptoncuts.cfg_min_pair_dca3d, dileptoncuts.cfg_max_pair_dca3d); // in sigma
    fDileptonCut.ApplyPhiV(dileptoncuts.cfg_apply_phiv);
    fDileptonCut.ApplyPrefilter(dileptoncuts.cfg_apply_pf);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    // for track
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_max_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMeanClusterSizeITSob(0, 16);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaMuRange(dileptoncuts.cfg_min_TPCNsigmaMu, dileptoncuts.cfg_max_TPCNsigmaMu);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTPCNsigmaKaRange(dileptoncuts.cfg_min_TPCNsigmaKa, dileptoncuts.cfg_max_TPCNsigmaKa);
    fDileptonCut.SetTPCNsigmaPrRange(dileptoncuts.cfg_min_TPCNsigmaPr, dileptoncuts.cfg_max_TPCNsigmaPr);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);

    if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      o2::ml::OnnxModel* eid_bdt = new o2::ml::OnnxModel();
      if (dileptoncuts.loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        std::map<std::string, std::string> metadata;
        bool retrieveSuccessGamma = ccdbApi.retrieveBlob(dileptoncuts.BDTPathCCDB.value, ".", metadata, dileptoncuts.timestampCCDB.value, false, dileptoncuts.BDTLocalPathGamma.value);
        if (retrieveSuccessGamma) {
          eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Gamma model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      } else {
        eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
      }

      fDileptonCut.SetPIDModel(eid_bdt);
    } // end of PID ML
  }

  template <typename TCollision, typename TTrack1, typename TTrack2>
  bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2)
  {
    if (t1.trackId() == t2.trackId()) { // this is protection against pairing identical 2 tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
      return false;
    }

    if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) {
      if (!fDileptonCut.IsSelectedTrack<true>(t1, collision) || !fDileptonCut.IsSelectedTrack<true>(t2, collision)) {
        return false;
      }
    } else { // cut-based
      if (!fDileptonCut.IsSelectedTrack(t1) || !fDileptonCut.IsSelectedTrack(t2)) {
        return false;
      }
    }

    if (!fDileptonCut.IsSelectedPair(t1, t2, collision.bz())) {
      return false;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    if (abs(v12.Rapidity()) > maxY) {
      return false;
    }
    float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision.bz());

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Pair/same/uls/hMvsPt"), v12.M(), v12.Pt());
      fRegistry.fill(HIST("Pair/same/uls/hMvsPhiV"), phiv, v12.M());
    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
      fRegistry.fill(HIST("Pair/same/lspp/hMvsPt"), v12.M(), v12.Pt());
      fRegistry.fill(HIST("Pair/same/lspp/hMvsPhiV"), phiv, v12.M());
    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
      fRegistry.fill(HIST("Pair/same/lsmm/hMvsPt"), v12.M(), v12.Pt());
      fRegistry.fill(HIST("Pair/same/lsmm/hMvsPhiV"), phiv, v12.M());
    }

    if (t1.sign() > 0) {
      if (std::find(used_trackIds.begin(), used_trackIds.end(), t1.globalIndex()) == used_trackIds.end()) {
        used_trackIds.emplace_back(t1.globalIndex());
        fillTrackInfo(t1);
      }
    } else {
      if (std::find(used_trackIds.begin(), used_trackIds.end(), t1.globalIndex()) == used_trackIds.end()) {
        used_trackIds.emplace_back(t1.globalIndex());
        fillTrackInfo(t1);
      }
    }
    if (t2.sign() > 0) {
      if (std::find(used_trackIds.begin(), used_trackIds.end(), t2.globalIndex()) == used_trackIds.end()) {
        used_trackIds.emplace_back(t2.globalIndex());
        fillTrackInfo(t2);
      }
    } else {
      if (std::find(used_trackIds.begin(), used_trackIds.end(), t2.globalIndex()) == used_trackIds.end()) {
        used_trackIds.emplace_back(t2.globalIndex());
        fillTrackInfo(t2);
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
    fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
    fRegistry.fill(HIST("Track/hDCAxyzSigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
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
  Filter trackFilter = static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && o2::aod::track::tpcChi2NCl < static_cast<float>(dileptoncuts.cfg_max_chi2tpc) && o2::aod::track::itsChi2NCl < static_cast<float>(dileptoncuts.cfg_max_chi2its) && nabs(o2::aod::track::dcaXY) < static_cast<float>(dileptoncuts.cfg_max_dcaxy) && nabs(o2::aod::track::dcaZ) < static_cast<float>(dileptoncuts.cfg_max_dcaz);
  Filter pidFilter = (static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl)) && (o2::aod::pidtpc::tpcNSigmaPi < static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaPi) || static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaPi) < o2::aod::pidtpc::tpcNSigmaPi);
  using FilteredMyTracks = soa::Filtered<MyTracks>;
  Partition<FilteredMyTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  std::vector<int> used_trackIds;
  void processQC(FilteredMyCollisions const& collisions, FilteredMyTracks const& /*tracks*/)
  {
    for (auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, cfgDoFlow);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, cfgDoFlow);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillPairInfo(collision, pos, ele);
      }
      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
        fillPairInfo(collision, pos1, pos2);
      }
      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
        fillPairInfo(collision, ele1, ele2);
      }
    } // end of collision loop

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();

  } // end of process
  PROCESS_SWITCH(DalitzEEQC, processQC, "run Dalitz QC", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(DalitzEEQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzEEQC>(cfgc, TaskName{"dalitz-ee-qc"})};
}
