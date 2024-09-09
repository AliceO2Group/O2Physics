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
// This code runs loop over v0 photons for PCM QC.
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct PCMQC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  EMPhotonEventCut fEMEventCut;
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

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<bool> cfg_require_v0_on_wwire_ib{"cfg_require_v0_on_wwire_ib", false, "flag to select V0s on W wires ITSib"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.997, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
    Configurable<bool> cfg_require_v0_with_correct_xz{"cfg_require_v0_with_correct_xz", true, "flag to select V0s with correct xz"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
  } pcmcuts;

  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    addhistograms();
    DefineEMEventCut();
    DefinePCMCut();
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

    // v0 info
    fRegistry.add("V0/hPt", "pT;p_{T,#gamma} (GeV/c)", kTH1F, {{2000, 0.0f, 20}}, false);
    fRegistry.add("V0/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, 2 * M_PI}, {40, -1.0f, 1.0f}}, false);
    fRegistry.add("V0/hRadius", "V0Radius; radius in Z (cm);radius in XY (cm)", kTH2F, {{200, -100, 100}, {200, 0.0f, 100.0f}}, false);
    fRegistry.add("V0/hCosPA", "V0CosPA;cosine pointing angle", kTH1F, {{100, 0.9f, 1.0f}}, false);
    fRegistry.add("V0/hCosPA_Rxy", "cos PA vs. R_{xy};R_{xy} (cm);cosine pointing angle", kTH2F, {{200, 0.f, 100.f}, {100, 0.9f, 1.0f}}, false);
    fRegistry.add("V0/hPCA", "distance between 2 legs;PCA (cm)", kTH1F, {{500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/hPCA_Rxy", "distance between 2 legs vs. R_{xy};R_{xy} (cm);PCA (cm)", kTH2F, {{200, 0.f, 100.f}, {500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/hPCA_CosPA", "distance between 2 legs vs. cosPA;cosine of pointing angle;PCA (cm)", kTH2F, {{100, 0.9f, 1.f}, {500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}, false);
    fRegistry.add("V0/hAPplot", "AP plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}}, false);
    fRegistry.add("V0/hMassGamma", "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.0f, 0.1f}}, false);
    fRegistry.add("V0/hGammaRxy", "conversion point in XY;V_{x} (cm);V_{y} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, -100.0f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsM", "KF chi2 vs. m_{ee};m_{ee} (GeV/c^{2});KF chi2/NDF", kTH2F, {{100, 0.0f, 0.1f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsR", "KF chi2 vs. conversion point in XY;R_{xy} (cm);KF chi2/NDF", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsX", "KF chi2 vs. conversion point in X;X (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsY", "KF chi2 vs. conversion point in Y;Y (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hKFChi2vsZ", "KF chi2 vs. conversion point in Z;Z (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/hNgamma", "Number of #gamma candidates per collision", kTH1F, {{101, -0.5f, 100.5f}});

    // v0leg info
    fRegistry.add("V0Leg/hPt", "pT;p_{T,e} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("V0Leg/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -50, 50}}, false);
    fRegistry.add("V0Leg/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, 2 * M_PI}, {40, -1.0f, 1.0f}}, false);
    fRegistry.add("V0Leg/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -50.0f, 50.0f}, {200, -50.0f, 50.0f}}, false);
    fRegistry.add("V0Leg/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("V0Leg/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("V0Leg/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("V0Leg/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {160, 0, 16}}, false);
    fRegistry.add("V0Leg/hXY", "X vs. Y;X (cm);Y (cm)", kTH2F, {{100, 0, 100}, {80, -20, 20}}, false);
    fRegistry.add("V0Leg/hZX", "Z vs. X;Z (cm);X (cm)", kTH2F, {{200, -100, 100}, {100, 0, 100}}, false);
    fRegistry.add("V0Leg/hZY", "Z vs. Y;Z (cm);Y (cm)", kTH2F, {{200, -100, 100}, {80, -20, 20}}, false);
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
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

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, 1e10f);
    fV0PhotonCut.SetV0EtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetTrackPtRange(pcmcuts.cfg_min_pt_v0 * 0.4, 1e+10f);
    fV0PhotonCut.SetTrackEtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfg_disable_itsonly_track);

    if (pcmcuts.cfg_reject_v0_on_itsib) {
      fV0PhotonCut.SetNClustersITS(2, 4);
    } else {
      fV0PhotonCut.SetNClustersITS(0, 7);
    }
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetIsWithinBeamPipe(pcmcuts.cfg_require_v0_with_correct_xz);

    if (pcmcuts.cfg_require_v0_with_itstpc) {
      fV0PhotonCut.SetRequireITSTPC(true);
      fV0PhotonCut.SetMaxPCA(1.0);
      fV0PhotonCut.SetRxyRange(4, 40);
    }
    if (pcmcuts.cfg_require_v0_with_itsonly) {
      fV0PhotonCut.SetRequireITSonly(true);
      fV0PhotonCut.SetMaxPCA(1.0);
      fV0PhotonCut.SetRxyRange(4, 24);
    }
    if (pcmcuts.cfg_require_v0_with_tpconly) {
      fV0PhotonCut.SetRequireTPConly(true);
      fV0PhotonCut.SetMaxPCA(3.0);
      fV0PhotonCut.SetRxyRange(36, 90);
    }
    if (pcmcuts.cfg_require_v0_on_wwire_ib) {
      fV0PhotonCut.SetMaxPCA(0.3);
      fV0PhotonCut.SetOnWwireIB(true);
      fV0PhotonCut.SetOnWwireOB(false);
      fV0PhotonCut.SetRxyRange(7, 14);
    }
  }

  template <const int ev_id, typename TCollision>
  void fillEventInfo(TCollision const& collision, const float /*weight*/ = 1.f)
  {
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
    }
    if (collision.sel8()) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
    }
    if (abs(collision.posZ()) < 10.0) {
      fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
    }
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0MvsMultNTracksPV"), collision.centFT0M(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0MvsMultNTracksPV"), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV());
  }

  template <typename TV0>
  void fillV0Info(TV0 const& v0)
  {
    fRegistry.fill(HIST("V0/hPt"), v0.pt());
    fRegistry.fill(HIST("V0/hEtaPhi"), v0.phi(), v0.eta());
    fRegistry.fill(HIST("V0/hRadius"), v0.vz(), v0.v0radius());
    fRegistry.fill(HIST("V0/hCosPA"), v0.cospa());
    fRegistry.fill(HIST("V0/hCosPA_Rxy"), v0.v0radius(), v0.cospa());
    fRegistry.fill(HIST("V0/hPCA"), v0.pca());
    fRegistry.fill(HIST("V0/hPCA_CosPA"), v0.cospa(), v0.pca());
    fRegistry.fill(HIST("V0/hPCA_Rxy"), v0.v0radius(), v0.pca());
    fRegistry.fill(HIST("V0/hDCAxyz"), v0.dcaXYtopv(), v0.dcaZtopv());
    fRegistry.fill(HIST("V0/hAPplot"), v0.alpha(), v0.qtarm());
    fRegistry.fill(HIST("V0/hMassGamma"), v0.v0radius(), v0.mGamma());
    fRegistry.fill(HIST("V0/hGammaRxy"), v0.vx(), v0.vy());
    fRegistry.fill(HIST("V0/hKFChi2vsM"), v0.mGamma(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsR"), v0.v0radius(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsX"), v0.vx(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsY"), v0.vy(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/hKFChi2vsZ"), v0.vz(), v0.chiSquareNDF());
  }

  template <typename TLeg>
  void fillV0LegInfo(TLeg const& leg)
  {
    fRegistry.fill(HIST("V0Leg/hPt"), leg.pt());
    fRegistry.fill(HIST("V0Leg/hQoverPt"), leg.sign() / leg.pt());
    fRegistry.fill(HIST("V0Leg/hEtaPhi"), leg.phi(), leg.eta());
    fRegistry.fill(HIST("V0Leg/hDCAxyz"), leg.dcaXY(), leg.dcaZ());
    fRegistry.fill(HIST("V0Leg/hNclsITS"), leg.itsNCls());
    fRegistry.fill(HIST("V0Leg/hNclsTPC"), leg.tpcNClsFound());
    fRegistry.fill(HIST("V0Leg/hNcrTPC"), leg.tpcNClsCrossedRows());
    fRegistry.fill(HIST("V0Leg/hTPCNcr2Nf"), leg.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("V0Leg/hTPCNcls2Nf"), leg.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("V0Leg/hChi2TPC"), leg.tpcChi2NCl());
    fRegistry.fill(HIST("V0Leg/hChi2ITS"), leg.itsChi2NCl());
    fRegistry.fill(HIST("V0Leg/hITSClusterMap"), leg.itsClusterMap());
    if (leg.hasITS()) {
      fRegistry.fill(HIST("V0Leg/hMeanClusterSizeITS"), leg.p(), leg.meanClusterSizeITS() * std::cos(std::atan(leg.tgl())));
    }
    fRegistry.fill(HIST("V0Leg/hTPCdEdx"), leg.tpcInnerParam(), leg.tpcSignal());
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaPi"), leg.tpcInnerParam(), leg.tpcNSigmaPi());
    fRegistry.fill(HIST("V0Leg/hXY"), leg.x(), leg.y());
    fRegistry.fill(HIST("V0Leg/hZX"), leg.z(), leg.x());
    fRegistry.fill(HIST("V0Leg/hZY"), leg.z(), leg.y());
  }

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy = eventcuts.cfgOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgOccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processQC(FilteredMyCollisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const&)
  {
    for (auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
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

      int nv0 = 0;
      auto v0photons_coll = v0photons.sliceBy(perCollision, collision.globalIndex());
      for (auto& v0 : v0photons_coll) {
        auto pos = v0.posTrack_as<aod::V0Legs>();
        auto ele = v0.negTrack_as<aod::V0Legs>();

        if (!fV0PhotonCut.IsSelected<aod::V0Legs>(v0)) {
          continue;
        }
        fillV0Info(v0);
        for (auto& leg : {pos, ele}) {
          fillV0LegInfo(leg);
        }
        nv0++;
      } // end of v0 loop
      fRegistry.fill(HIST("V0/hNgamma"), nv0);
    } // end of collision loop
  } // end of process

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(PCMQC, processQC, "run PCM QC", true);
  PROCESS_SWITCH(PCMQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQC>(cfgc, TaskName{"pcm-qc"})};
}
