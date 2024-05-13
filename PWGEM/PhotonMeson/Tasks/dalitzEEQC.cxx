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
#include "TString.h"
#include "THashList.h"
#include "TDirectory.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMEventCut.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EMTrack.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;
using std::array;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsBz>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;

struct DalitzEEQC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};

  EMEventCut fEMEventCut;
  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
  Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
  Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
  Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
  Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
  Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
  Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};

  DalitzEECut fDielectronCut;
  Configurable<float> cfg_min_mee{"cfg_min_mee", 0.0, "min mee"};
  Configurable<float> cfg_max_mee{"cfg_max_mee", 1e+10, "max mee"};
  Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
  Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
  Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
  Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};

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
  Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron inclusion"};
  Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
  Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
  Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
  Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -3.0, "min. TPC n sigma for pion exclusion"};
  Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
  Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
  Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
  Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
  Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
  Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> cent_bin_edges;
  std::vector<float> zvtx_bin_edges;

  void init(InitContext& context)
  {
    zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
    zvtx_bin_edges.erase(zvtx_bin_edges.begin());

    cent_bin_edges = std::vector<float>(ConfCentBins.value.begin(), ConfCentBins.value.end());
    cent_bin_edges.erase(cent_bin_edges.begin());

    DefineEMEventCut();
    DefineDielectronCut();
    addhistograms(); // please call this after DefinCuts();
  }

  ~DalitzEEQC()
  {
    map_mix_bins.clear();
    map_posTracks_to_collision.clear();
    map_negTracks_to_collision.clear();

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();
  }

  void addhistograms()
  {
    // event info
    auto hCollisionCounter = fRegistry.add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{10, 0.5, 10.5}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(10, "accepted");

    fRegistry.add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{300, 0, 6000}, {300, 0, 6000}}, false);
    fRegistry.add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/before/hCentFT0MvsMultNTracksPV", "hCentFT0MvsMultNTracksPV;centrality FT0M (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
    fRegistry.add("Event/before/hMultFT0MvsMultNTracksPV", "hMultFT0MvsMultNTracksPV;mult. FT0M;N_{track} to PV", kTH2F, {{600, 0, 6000}, {600, 0, 6000}}, false);
    fRegistry.addClone("Event/before/", "Event/after/");

    // pair info

    const int nm = 145;
    double mee[nm] = {0.f};
    for (int i = 0; i < 110; i++) {
      mee[i] = 0.01 * (i - 0) + 0.0; // every 0.01 GeV/c2 up to 1.1 GeV/c2
    }
    for (int i = 110; i < 126; i++) {
      mee[i] = 0.1 * (i - 110) + 1.1; // every 0.1 GeV/c2 from 1.1 to 2.7 GeV/c2
    }
    for (int i = 126; i < 136; i++) {
      mee[i] = 0.05 * (i - 126) + 2.7; // every 0.05 GeV/c2 from 2.7 to 3.2 GeV/c2
    }
    for (int i = 136; i < nm; i++) {
      mee[i] = 0.1 * (i - 136) + 3.2; // every 0.1 GeV/c2 from 3.2 to 4.0 GeV/c2
    }

    const int npt = 61;
    double pt[npt] = {0.f};
    for (int i = 0; i < 50; i++) {
      pt[i] = 0.1 * i;
    }
    for (int i = 50; i < npt; i++) {
      pt[i] = 0.5 * (i - 50) + 5.0;
    }

    const int ndca = 27;
    double dca[ndca] = {0.f};
    for (int i = 0; i < 20; i++) {
      dca[i] = 0.1 * i;
    }
    for (int i = 20; i < ndca; i++) {
      dca[i] = 0.5 * (i - 20) + 2.0;
    }

    auto hs_uls = fRegistry.add<THnSparse>("Pair/same/hs_uls", "hs pair;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);DCA_{ee}^{3D} (#sigma);", kTHnSparseF, {{nm - 1, 0, 4.0}, {npt - 1, 0, 10}, {ndca - 1, 0, 5}}, true);
    hs_uls->SetBinEdges(0, mee);
    hs_uls->SetBinEdges(1, pt);
    hs_uls->SetBinEdges(2, dca);
    fRegistry.addClone("Pair/same/hs_uls", "Pair/same/hs_lspp");
    fRegistry.addClone("Pair/same/hs_uls", "Pair/same/hs_lsmm");

    fRegistry.add("Pair/same/hMvsPhiV_uls", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, false);
    fRegistry.addClone("Pair/same/hMvsPhiV_uls", "Pair/same/hMvsPhiV_lspp");
    fRegistry.addClone("Pair/same/hMvsPhiV_uls", "Pair/same/hMvsPhiV_lsmm");
    fRegistry.addClone("Pair/same/", "Pair/mix/");

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
    fEMEventCut.SetRequireSel8(cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-cfgZvtxMax, +cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(cfgRequireGoodZvtxFT0vsPV);
  }

  void DefineDielectronCut()
  {
    fDielectronCut = DalitzEECut("fDielectronCut", "fDielectronCut");

    // for pair
    fDielectronCut.SetMeeRange(cfg_min_mee, cfg_max_mee);
    fDielectronCut.SetMaxPhivPairMeeDep([](float mee) { return (mee - -0.028) / 0.0185; });
    fDielectronCut.SetPairDCARange(cfg_min_pair_dca3d, cfg_max_pair_dca3d); // in sigma
    fDielectronCut.ApplyPhiV(cfg_apply_phiv);
    fDielectronCut.ApplyPrefilter(cfg_apply_pf);

    // for track
    fDielectronCut.SetTrackPtRange(cfg_min_pt_track, 1e+10f);
    fDielectronCut.SetTrackEtaRange(-cfg_max_eta_track, +cfg_max_eta_track);
    fDielectronCut.SetMinNClustersTPC(cfg_min_ncluster_tpc);
    fDielectronCut.SetMinNCrossedRowsTPC(cfg_min_ncrossedrows);
    fDielectronCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDielectronCut.SetChi2PerClusterTPC(0.0, cfg_max_chi2tpc);
    fDielectronCut.SetChi2PerClusterITS(0.0, cfg_max_chi2its);
    fDielectronCut.SetNClustersITS(cfg_min_ncluster_its, 7);
    fDielectronCut.SetMeanClusterSizeITSob(0, 16);
    fDielectronCut.SetMaxDcaXY(cfg_max_dcaxy);
    fDielectronCut.SetMaxDcaZ(cfg_max_dcaz);

    // for eID
    fDielectronCut.SetPIDScheme(cfg_pid_scheme);
    fDielectronCut.SetTPCNsigmaElRange(cfg_min_TPCNsigmaEl, cfg_max_TPCNsigmaEl);
    fDielectronCut.SetTPCNsigmaMuRange(cfg_min_TPCNsigmaMu, cfg_max_TPCNsigmaMu);
    fDielectronCut.SetTPCNsigmaPiRange(cfg_min_TPCNsigmaPi, cfg_max_TPCNsigmaPi);
    fDielectronCut.SetTPCNsigmaKaRange(cfg_min_TPCNsigmaKa, cfg_max_TPCNsigmaKa);
    fDielectronCut.SetTPCNsigmaPrRange(cfg_min_TPCNsigmaPr, cfg_max_TPCNsigmaPr);
    fDielectronCut.SetTOFNsigmaElRange(cfg_min_TOFNsigmaEl, cfg_max_TOFNsigmaEl);
  }

  template <const int ev_id, typename TCollision>
  void fillEventInfo(TCollision const& collision, const float weight = 1.f)
  {
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
    }
    if (collision.sel8()) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
    }
    if (abs(collision.posZ()) < 10.0) {
      fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
    }
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hZvtx"), collision.posZ());

    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hCentFT0MvsMultNTracksPV"), collision.centFT0M(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_cut_types[ev_id]) + HIST("hMultFT0MvsMultNTracksPV"), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV());
  }

  template <int ev_id, typename TCollision, typename TTrack1, typename TTrack2>
  bool fillPairInfo(TCollision const& collision, TTrack1 const& t1, TTrack2 const& t2, const float weight = 1.f)
  {
    if (t1.trackId() == t2.trackId()) { // this is protection against pairing identical 2 tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
      return false;
    }

    if constexpr (ev_id == 0) {
      if (!fDielectronCut.IsSelectedTrack(t1) || !fDielectronCut.IsSelectedTrack(t2)) {
        return false;
      }
    }

    if (!fDielectronCut.IsSelectedPair(t1, t2, collision.bz())) {
      return false;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dca_t1_3d = t1.dca3DinSigma();
    float dca_t2_3d = t2.dca3DinSigma();
    float dca_ee_3d = std::sqrt((dca_t1_3d * dca_t1_3d + dca_t2_3d * dca_t2_3d) / 2.);
    float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision.bz());

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_uls"), v12.M(), v12.Pt(), dca_ee_3d);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hMvsPhiV_uls"), phiv, v12.M());
    } else if (t1.sign() > 0 && t2.sign() > 0) { // LS++
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_lspp"), v12.M(), v12.Pt(), dca_ee_3d);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hMvsPhiV_lspp"), phiv, v12.M());
    } else if (t1.sign() < 0 && t2.sign() < 0) { // LS--
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_lsmm"), v12.M(), v12.Pt(), dca_ee_3d);
      fRegistry.fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hMvsPhiV_lsmm"), phiv, v12.M());
    }

    // store tracks for event mixing without double counting
    if constexpr (ev_id == 0) {
      if (t1.sign() > 0) {
        if (std::find(used_trackIds.begin(), used_trackIds.end(), std::make_pair(ndf, t1.trackId())) == used_trackIds.end()) {
          map_posTracks_to_collision[std::make_pair(ndf, collision.globalIndex())].emplace_back(EMTrack(collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), t1.sign(), t1.dca3DinSigma()));
          used_trackIds.emplace_back(std::make_pair(ndf, t1.trackId()));
        }
      } else {
        if (std::find(used_trackIds.begin(), used_trackIds.end(), std::make_pair(ndf, t1.trackId())) == used_trackIds.end()) {
          map_negTracks_to_collision[std::make_pair(ndf, collision.globalIndex())].emplace_back(EMTrack(collision.globalIndex(), t1.trackId(), t1.pt(), t1.eta(), t1.phi(), t1.sign(), t1.dca3DinSigma()));
          used_trackIds.emplace_back(std::make_pair(ndf, t1.trackId()));
        }
      }
      if (t2.sign() > 0) {
        if (std::find(used_trackIds.begin(), used_trackIds.end(), std::make_pair(ndf, t2.trackId())) == used_trackIds.end()) {
          map_posTracks_to_collision[std::make_pair(ndf, collision.globalIndex())].emplace_back(EMTrack(collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), t2.sign(), t2.dca3DinSigma()));
          used_trackIds.emplace_back(std::make_pair(ndf, t2.trackId()));
        }
      } else {
        if (std::find(used_trackIds.begin(), used_trackIds.end(), std::make_pair(ndf, t2.trackId())) == used_trackIds.end()) {
          map_negTracks_to_collision[std::make_pair(ndf, collision.globalIndex())].emplace_back(EMTrack(collision.globalIndex(), t2.trackId(), t2.pt(), t2.eta(), t2.phi(), t2.sign(), t2.dca3DinSigma()));
          used_trackIds.emplace_back(std::make_pair(ndf, t2.trackId()));
        }
      }
    }
    return true;
  }

  template <typename TTrack>
  void fillTrackInfo(TTrack const& track, const float weight = 1.f)
  {
    if (!fDielectronCut.IsSelectedTrack(track)) {
      return;
    }

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

  std::map<std::tuple<int, int, int>, std::vector<std::pair<int, int64_t>>> map_mix_bins; // zvtx, centrality, event plane bins -> pair<df index, vector of collision>
  std::map<std::pair<int, int64_t>, std::vector<EMTrack>> map_posTracks_to_collision;     // pair<df index, collisionId> -> vector of track
  std::map<std::pair<int, int64_t>, std::vector<EMTrack>> map_negTracks_to_collision;     // pair<df index, collisionId> -> vector of track

  SliceCache cache;
  Preslice<MyTracks> perCollision_track = aod::emprimaryelectron::emeventId;
  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  Filter trackFilter = cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < cfg_max_eta_track && (cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < cfg_max_TPCNsigmaEl);
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<MyFilteredTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  std::vector<std::pair<int, int>> used_trackIds;
  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {VARIABLE_WIDTH, 0.0f, M_PI / 4, M_PI / 2, M_PI}, "Mixing bins - event plane angle"};

  int ndf = 0;
  void processQC(MyCollisions const& collisions, MyFilteredTracks const& tracks)
  {
    for (auto& collision : grouped_collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      fillEventInfo<0>(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      fillEventInfo<1>(collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

      for (auto& pos : posTracks_per_coll) {
        fillTrackInfo(pos);
      }

      for (auto& ele : negTracks_per_coll) {
        fillTrackInfo(ele);
      }

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

      // make a vector of selected tracks in this collision.
      auto selected_posTracks_in_this_event = map_posTracks_to_collision[std::make_pair(ndf, collision.globalIndex())];
      auto selected_negTracks_in_this_event = map_negTracks_to_collision[std::make_pair(ndf, collision.globalIndex())];
      // LOGF(info, "N selected tracks in current event (%d, %d)", selected_posTracks_in_this_event.size(), selected_negTracks_in_this_event.size());

      auto collisionIds_in_mixing_pool = map_mix_bins[std::make_tuple(zbin, centbin, 0)];
      for (auto& mix_dfId_collisionId : collisionIds_in_mixing_pool) {
        int mix_dfId = mix_dfId_collisionId.first;
        int64_t mix_collisionId = mix_dfId_collisionId.second;

        if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
          continue;
        }

        auto posTracks_from_event_pool = map_posTracks_to_collision[mix_dfId_collisionId];
        auto negTracks_from_event_pool = map_negTracks_to_collision[mix_dfId_collisionId];
        // LOGF(info, "Do event mixing: current event (%d, %d) | event pool (%d, %d), npos = %d , nele = %d", ndf, collision.globalIndex(), mix_dfId, mix_collisionId, posTracks_from_event_pool.size(), negTracks_from_event_pool.size());

        for (auto& pos : selected_posTracks_in_this_event) { // ULS mix
          for (auto& ele : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos, ele);
          }
        }

        for (auto& ele : selected_negTracks_in_this_event) { // ULS mix
          for (auto& pos : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, ele, pos);
          }
        }

        for (auto& pos1 : selected_posTracks_in_this_event) { // LS++ mix
          for (auto& pos2 : posTracks_from_event_pool) {
            fillPairInfo<1>(collision, pos1, pos2);
          }
        }

        for (auto& ele1 : selected_negTracks_in_this_event) { // LS-- mix
          for (auto& ele2 : negTracks_from_event_pool) {
            fillPairInfo<1>(collision, ele1, ele2);
          }
        }
      } // end of loop over mixed event pool

      if (nuls > 0 || nlspp > 0 || nlsmm > 0) { // prepare fill bins
        if (static_cast<int>(map_mix_bins[std::make_tuple(zbin, centbin, 0)].size()) >= ndepth) {
          // LOGF(info, "nev = %d , erase elements.",  map_mix_bins[std::make_tuple(zbin, centbin, 0)].size() );
          map_mix_bins[std::make_tuple(zbin, centbin, 0)].erase(map_mix_bins[std::make_tuple(zbin, centbin, 0)].begin());
          // LOGF(info, "after erase : nev = %d",  map_mix_bins[std::make_tuple(zbin, centbin, 0)].size() );

          map_posTracks_to_collision[map_mix_bins[std::make_tuple(zbin, centbin, 0)][0]].clear();
          map_posTracks_to_collision[map_mix_bins[std::make_tuple(zbin, centbin, 0)][0]].shrink_to_fit();
          map_negTracks_to_collision[map_mix_bins[std::make_tuple(zbin, centbin, 0)][0]].clear();
          map_negTracks_to_collision[map_mix_bins[std::make_tuple(zbin, centbin, 0)][0]].shrink_to_fit();
        }
        // LOGF(info, "npos = %d , nele = %d after remove",  map_posTracks_to_collision[map_mix_bins[std::make_tuple(zbin, centbin, 0)][0]].size(),map_negTracks_to_collision[map_mix_bins[std::make_tuple(zbin, centbin, 0)][0]].size() );
        map_mix_bins[std::make_tuple(zbin, centbin, 0)].emplace_back(std::make_pair(ndf, collision.globalIndex()));
        // LOGF(info, "after emplace_back : nev = %d",  map_mix_bins[std::make_tuple(zbin, centbin, 0)].size() );
      }
    } // end of collision loop

    ndf++;
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
