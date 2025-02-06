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

#include <string>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::BinnedGenPts>;
using MyMCCollision = MyMCCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0PhotonsKFCov, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

struct PCMQCMC {

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<float> maxRgen{"maxRgen", 90.f, "maximum radius for generated particles"};
  Configurable<float> margin_z_mc{"margin_z_mc", 7.0, "margin for z cut in cm for MC"};

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

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<bool> cfg_require_v0_on_wwire_ib{"cfg_require_v0_on_wwire_ib", false, "flag to select V0s on W wires ITSib"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", +0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.997, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
    Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfg_require_v0_with_correct_xz{"cfg_require_v0_with_correct_xz", true, "flag to select V0s with correct xz"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
  } pcmcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  static constexpr std::string_view mcphoton_types[3] = {"primary/", "fromWD/", "fromHS/"};

  void init(InitContext&)
  {
    DefineEMEventCut();
    DefinePCMCut();
    addhistograms();
  }

  void addhistograms()
  {
    std::vector<double> ptbins;
    for (int i = 0; i < 2; i++) {
      ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.05 GeV/c, every 0.05 GeV/c
    }
    for (int i = 2; i < 51; i++) {
      ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 4.9 GeV/c, every 0.1 GeV/c
    }
    for (int i = 51; i < 61; i++) {
      ptbins.emplace_back(0.5 * (i - 51) + 5.0); // from 5 to 9.5 GeV/c, every 0.5 GeV/c
    }
    for (int i = 61; i < 72; i++) {
      ptbins.emplace_back(1.0 * (i - 61) + 10.0); // from 10 to 20 GeV/c, every 1 GeV/c
    }
    const AxisSpec axis_pt{ptbins, "p_{T,#gamma} (GeV/c)"};
    const AxisSpec axis_rapidity{{0.0, +0.8, +0.9}, "rapidity |y_{#gamma}|"};

    if (doprocessGen) {
      fRegistry.add("Generated/hPt", "pT;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);
      fRegistry.add("Generated/hPtY", "Generated info", kTH2F, {axis_pt, axis_rapidity}, true);
      fRegistry.add("Generated/hPt_ConvertedPhoton", "converted photon pT;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);
      fRegistry.add("Generated/hY_ConvertedPhoton", "converted photon y;rapidity y", kTH1F, {{40, -2.0f, 2.0f}}, true);
      fRegistry.add("Generated/hPhi_ConvertedPhoton", "converted photon #varphi;#varphi (rad.)", kTH1F, {{180, 0, 2 * M_PI}}, true);
      fRegistry.add("Generated/hPhotonRxy", "conversion point in XY MC;V_{x} (cm);V_{y} (cm)", kTH2F, {{800, -100.0f, 100.0f}, {800, -100.0f, 100.0f}}, true);
      fRegistry.add("Generated/hPhotonRZ", "conversion point in RZ MC;V_{z} (cm);R_{xy} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, 0.f, 100.0f}}, true);
      fRegistry.add("Generated/hPhotonPhivsRxy", "conversion point of #varphi vs. R_{xy} MC;#varphi (rad.);R_{xy} (cm);N_{e}", kTH2F, {{360, 0.0f, 2 * M_PI}, {400, 0, 100}}, true);
    }

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
    fRegistry.add("V0/primary/hPt", "pT;p_{T,#gamma} (GeV/c)", kTH1F, {{2000, 0.0f, 20}}, false);
    fRegistry.add("V0/primary/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, 2 * M_PI}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("V0/primary/hRadius", "V0Radius; radius in Z (cm);radius in XY (cm)", kTH2F, {{200, -100, 100}, {200, 0.0f, 100.0f}}, false);
    fRegistry.add("V0/primary/hCosPA", "V0CosPA;cosine pointing angle", kTH1F, {{100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/primary/hCosPA_Rxy", "cos PA vs. R_{xy};R_{xy} (cm);cosine pointing angle", kTH2F, {{200, 0.f, 100.f}, {100, 0.99f, 1.0f}}, false);
    fRegistry.add("V0/primary/hPCA", "distance between 2 legs;PCA (cm)", kTH1F, {{500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/primary/hPCA_Rxy", "distance between 2 legs vs. R_{xy};R_{xy} (cm);PCA (cm)", kTH2F, {{200, 0.f, 100.f}, {500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/primary/hPCA_CosPA", "distance between 2 legs vs. cosPA;cosine of pointing angle;PCA (cm)", kTH2F, {{100, 0.99f, 1.f}, {500, 0.0f, 5.0f}}, false);
    fRegistry.add("V0/primary/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}, false);
    fRegistry.add("V0/primary/hAPplot", "AP plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}}, false);
    fRegistry.add("V0/primary/hMassGamma", "hMassGamma;R_{xy} (cm);m_{ee} (GeV/c^{2})", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.0f, 0.1f}}, false);
    fRegistry.add("V0/primary/hGammaRxy", "conversion point in XY;V_{x} (cm);V_{y} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, -100.0f, 100.0f}}, false);
    fRegistry.add("V0/primary/hKFChi2vsM", "KF chi2 vs. m_{ee};m_{ee} (GeV/c^{2});KF chi2/NDF", kTH2F, {{100, 0.0f, 0.1f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/primary/hKFChi2vsR", "KF chi2 vs. conversion point in XY;R_{xy} (cm);KF chi2/NDF", kTH2F, {{200, 0.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/primary/hKFChi2vsX", "KF chi2 vs. conversion point in X;X (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/primary/hKFChi2vsY", "KF chi2 vs. conversion point in Y;Y (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/primary/hKFChi2vsZ", "KF chi2 vs. conversion point in Z;Z (cm);KF chi2/NDF", kTH2F, {{200, -100.0f, 100.0f}, {100, 0.f, 100.0f}}, false);
    fRegistry.add("V0/primary/hPResolution", "p resolution;p_{#gamma} (GeV/c);#Deltap/p", kTH2F, {{1000, 0.0f, 10}, {100, 0, 0.1}}, false);
    fRegistry.add("V0/primary/hPtResolution", "p_{T} resolution;p_{#gamma} (GeV/c);#Deltap_{T}/p_{T}", kTH2F, {{1000, 0.0f, 10}, {100, 0, 0.1}}, false);
    fRegistry.add("V0/primary/hEtaResolution", "#eta resolution;p_{#gamma} (GeV/c);#Delta#eta", kTH2F, {{1000, 0.0f, 10}, {100, 0, 0.01}}, false);
    fRegistry.add("V0/primary/hThetaResolution", "#theta resolution;p_{#gamma} (GeV/c);#Delta#theta (rad.)", kTH2F, {{1000, 0.0f, 10}, {100, 0, 0.01}}, false);
    fRegistry.add("V0/primary/hPhiResolution", "#varphi resolution;p_{#gamma} (GeV/c);#Delta#varphi (rad.)", kTH2F, {{1000, 0.0f, 10}, {100, 0, 0.01}}, false);
    fRegistry.add("V0/primary/hNgamma", "Number of true #gamma per collision", kTH1F, {{101, -0.5f, 100.5f}});
    fRegistry.add("V0/primary/hConvPoint_diffX", "conversion point diff X MC;X_{MC} (cm);X_{rec} - X_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
    fRegistry.add("V0/primary/hConvPoint_diffY", "conversion point diff Y MC;Y_{MC} (cm);Y_{rec} - Y_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
    fRegistry.add("V0/primary/hConvPoint_diffZ", "conversion point diff Z MC;Z_{MC} (cm);Z_{rec} - Z_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
    fRegistry.add("V0/primary/hPtGen_DeltaPtOverPtGen", "photon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/primary/hPtGen_DeltaEta", "photon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/primary/hPtGen_DeltaPhi", "photon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/primary/hXY_Photon_MC", "X vs. Y of true photon conversion point.;X (cm);Y (cm)", kTH2F, {{400, -100.0f, +100}, {400, -100, +100}}, true);
    fRegistry.add("V0/primary/hRZ_Photon_MC", "R vs. Z of true photon conversion point;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100.0f, +100}, {200, 0, 100}}, true);
    fRegistry.addClone("V0/primary/", "V0/fromWD/"); // from weak decay
    fRegistry.addClone("V0/primary/", "V0/fromHS/"); // from hadronic shower in detector materials

    // v0leg info
    fRegistry.add("V0Leg/primary/hPt", "pT;p_{T,e} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("V0Leg/primary/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -50, 50}}, false);
    fRegistry.add("V0Leg/primary/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, 2 * M_PI}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("V0Leg/primary/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -50.0f, 50.0f}, {200, -50.0f, 50.0f}}, false);
    fRegistry.add("V0Leg/primary/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/primary/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("V0Leg/primary/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/primary/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("V0Leg/primary/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/primary/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("V0Leg/primary/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/primary/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0Leg/primary/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    fRegistry.add("V0Leg/primary/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("V0Leg/primary/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0Leg/primary/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("V0Leg/primary/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {160, 0, 16}}, false);
    fRegistry.add("V0Leg/primary/hXY", "X vs. Y;X (cm);Y (cm)", kTH2F, {{100, 0, 100}, {40, -20, 20}}, false);
    fRegistry.add("V0Leg/primary/hZX", "Z vs. X;Z (cm);X (cm)", kTH2F, {{200, -100, 100}, {100, 0, 100}}, false);
    fRegistry.add("V0Leg/primary/hZY", "Z vs. Y;Z (cm);Y (cm)", kTH2F, {{200, -100, 100}, {40, -20, 20}}, false);
    fRegistry.add("V0Leg/primary/hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0Leg/primary/hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0Leg/primary/hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.addClone("V0Leg/primary/", "V0Leg/fromWD/"); // from weak decay
    fRegistry.addClone("V0Leg/primary/", "V0Leg/fromHS/"); // from hadronic shower in detector materials
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
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, 1e10f);
    fV0PhotonCut.SetV0EtaRange(pcmcuts.cfg_min_eta_v0, pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfg_max_chi2kf);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetTrackPtRange(pcmcuts.cfg_min_pt_v0 * 0.4, 1e+10f);
    fV0PhotonCut.SetTrackEtaRange(pcmcuts.cfg_min_eta_v0, pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfg_max_frac_shared_clusters_tpc);
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
      fV0PhotonCut.SetRxyRange(4, 40);
    }
    if (pcmcuts.cfg_require_v0_with_itsonly) {
      fV0PhotonCut.SetRequireITSonly(true);
      fV0PhotonCut.SetRxyRange(4, 24);
    }
    if (pcmcuts.cfg_require_v0_with_tpconly) {
      fV0PhotonCut.SetRequireTPConly(true);
      fV0PhotonCut.SetRxyRange(32, 90);
    }
    if (pcmcuts.cfg_require_v0_on_wwire_ib) {
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

  template <int mctype, typename TV0, typename TMCV0, typename TMCLeg>
  void fillV0Info(TV0 const& v0, TMCV0 const& mcphoton, TMCLeg const& mcleg)
  {
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPt"), v0.pt());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hEtaPhi"), v0.phi(), v0.eta());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRadius"), v0.vz(), v0.v0radius());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hCosPA"), v0.cospa());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hCosPA_Rxy"), v0.v0radius(), v0.cospa());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPCA"), v0.pca());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPCA_CosPA"), v0.cospa(), v0.pca());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPCA_Rxy"), v0.v0radius(), v0.pca());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hDCAxyz"), v0.dcaXYtopv(), v0.dcaZtopv());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hAPplot"), v0.alpha(), v0.qtarm());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hMassGamma"), v0.v0radius(), v0.mGamma());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hGammaRxy"), v0.vx(), v0.vy());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsM"), v0.mGamma(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsR"), v0.v0radius(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsX"), v0.vx(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsY"), v0.vy(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hKFChi2vsZ"), v0.vz(), v0.chiSquareNDF());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPResolution"), v0.p(), getPResolution(v0) / v0.p());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtResolution"), v0.p(), getPtResolution(v0) / v0.pt());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hEtaResolution"), v0.p(), getEtaResolution(v0));
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hThetaResolution"), v0.p(), getThetaResolution(v0));
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPhiResolution"), v0.p(), getPhiResolution(v0));
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPtOverPtGen"), mcphoton.pt(), (v0.pt() - mcphoton.pt()) / mcphoton.pt());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaEta"), mcphoton.pt(), v0.eta() - mcphoton.eta());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPhi"), mcphoton.pt(), v0.phi() - mcphoton.phi());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hConvPoint_diffX"), mcleg.vx(), v0.vx() - mcleg.vx());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hConvPoint_diffY"), mcleg.vy(), v0.vy() - mcleg.vy());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hConvPoint_diffZ"), mcleg.vz(), v0.vz() - mcleg.vz());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hXY_Photon_MC"), mcleg.vx(), mcleg.vy());
    fRegistry.fill(HIST("V0/") + HIST(mcphoton_types[mctype]) + HIST("hRZ_Photon_MC"), mcleg.vz(), std::sqrt(std::pow(mcleg.vx(), 2) + std::pow(mcleg.vy(), 2)));
  }

  template <int mctype, typename TLeg>
  void fillV0LegInfo(TLeg const& leg)
  {
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPt"), leg.pt());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hQoverPt"), leg.sign() / leg.pt());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hEtaPhi"), leg.phi(), leg.eta());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hDCAxyz"), leg.dcaXY(), leg.dcaZ());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hNclsITS"), leg.itsNCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hNclsTPC"), leg.tpcNClsFound());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hNcrTPC"), leg.tpcNClsCrossedRows());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNcr2Nf"), leg.tpcCrossedRowsOverFindableCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNcls2Nf"), leg.tpcFoundOverFindableCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNclsShared"), leg.pt(), leg.tpcFractionSharedCls());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hChi2TPC"), leg.tpcChi2NCl());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hChi2ITS"), leg.itsChi2NCl());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hITSClusterMap"), leg.itsClusterMap());
    if (leg.hasITS()) {
      fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hMeanClusterSizeITS"), leg.p(), leg.meanClusterSizeITS() * std::cos(std::atan(leg.tgl())));
    }
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCdEdx"), leg.tpcInnerParam(), leg.tpcSignal());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hTPCNsigmaPi"), leg.tpcInnerParam(), leg.tpcNSigmaPi());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hXY"), leg.x(), leg.y());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hZX"), leg.z(), leg.x());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hZY"), leg.z(), leg.y());
    auto mcleg = leg.template emmcparticle_as<aod::EMMCParticles>();
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPtOverPtGen"), mcleg.pt(), (leg.pt() - mcleg.pt()) / mcleg.pt());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaEta"), mcleg.pt(), leg.eta() - mcleg.eta());
    fRegistry.fill(HIST("V0Leg/") + HIST(mcphoton_types[mctype]) + HIST("hPtGen_DeltaPhi"), mcleg.pt(), leg.phi() - mcleg.phi());
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_occupancy_track = eventcuts.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventcuts.cfgTrackOccupancyMax;
  Filter collisionFilter_occupancy_ft0c = eventcuts.cfgFT0COccupancyMin < o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventcuts.cfgFT0COccupancyMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  void processQCMC(FilteredMyCollisions const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const&, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const&)
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

      auto V0Photons_coll = v0photons.sliceBy(perCollision, collision.globalIndex());
      int ng_primary = 0, ng_wd = 0, ng_hs = 0;
      for (auto& v0 : V0Photons_coll) {
        auto pos = v0.posTrack_as<MyMCV0Legs>();
        auto ele = v0.negTrack_as<MyMCV0Legs>();
        auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
        auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

        // LOGF(info, "posmc.isPhysicalPrimary() = %d, posmc.producedByGenerator() = %d, elemc.isPhysicalPrimary() = %d, elemc.producedByGenerator() = %d", posmc.isPhysicalPrimary(), posmc.producedByGenerator(), elemc.isPhysicalPrimary(), elemc.producedByGenerator());

        if (!fV0PhotonCut.IsSelected<MyMCV0Legs>(v0)) {
          continue;
        }
        int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
        if (photonid < 0) {
          continue;
        }
        auto mcphoton = mcparticles.iteratorAt(photonid);

        if (mcphoton.isPhysicalPrimary() || mcphoton.producedByGenerator()) {
          fillV0Info<0>(v0, mcphoton, elemc);
          for (auto& leg : {pos, ele}) {
            fillV0LegInfo<0>(leg);
          }
          ng_primary++;
        } else if (IsFromWD(mcphoton.template emmcevent_as<aod::EMMCEvents>(), mcphoton, mcparticles) > 0) {
          fillV0Info<1>(v0, mcphoton, elemc);
          for (auto& leg : {pos, ele}) {
            fillV0LegInfo<1>(leg);
          }
          ng_wd++;
        } else {
          fillV0Info<2>(v0, mcphoton, elemc);
          for (auto& leg : {pos, ele}) {
            fillV0LegInfo<2>(leg);
          }
          ng_hs++;
          // LOGF(info, "mcphoton.vx() = %f, mcphoton.vy() = %f, mcphoton.vz() = %f, mother_pdg = %d", mcphoton.vx(), mcphoton.vy(), mcphoton.vz(), mother_pdg);
        }
      } // end of v0 loop
      fRegistry.fill(HIST("V0/primary/hNgamma"), ng_primary);
      fRegistry.fill(HIST("V0/fromWD/hNgamma"), ng_wd);
      fRegistry.fill(HIST("V0/fromHS/hNgamma"), ng_hs);
    } // end of collision loop
  } // end of process

  template <typename TBinnedData>
  void fillBinnedData(TBinnedData const& binned_data, const float weight = 1.f)
  {
    int xbin = 0, ybin = 0, zbin = 0;
    auto hPtY = fRegistry.get<TH2>(HIST("Generated/hPtY")); // 2D
    auto hPt = fRegistry.get<TH1>(HIST("Generated/hPt"));   // 1D

    for (int ibin = 0; ibin < hPtY->GetNcells(); ibin++) {
      int nentry = binned_data[ibin];
      hPtY->GetBinXYZ(ibin, xbin, ybin, zbin);
      float pt = hPtY->GetXaxis()->GetBinCenter(xbin);
      float y = hPtY->GetYaxis()->GetBinCenter(ybin);
      if (y > pcmcuts.cfg_max_eta_v0) {
        continue;
      }

      for (int j = 0; j < nentry; j++) {
        hPtY->Fill(pt, y, weight);
        hPt->Fill(pt, weight);
      }
    }
  }

  // aod::EMMCParticles contain primary photon conversions. // aod::BinnedGenPts contain primary photons
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(FilteredMyCollisions const& collisions, MyMCCollisions const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto mccollision = collision.template emmcevent_as<MyMCCollisions>();
      auto binned_data_gamma_gen = mccollision.generatedGamma();
      fillBinnedData(binned_data_gamma_gen, 1.f);
      // LOGF(info, "mccollision.globalIndex() = %d", mccollision.globalIndex());

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {
        if (abs(mctrack.y()) > pcmcuts.cfg_max_eta_v0) {
          continue;
        }

        if (abs(mctrack.pdgCode()) == 22 && (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          fRegistry.fill(HIST("Generated/hPt_ConvertedPhoton"), mctrack.pt());
          fRegistry.fill(HIST("Generated/hY_ConvertedPhoton"), mctrack.y());
          fRegistry.fill(HIST("Generated/hPhi_ConvertedPhoton"), mctrack.phi());
          auto daughter = mcparticles.iteratorAt(mctrack.daughtersIds()[0]); // choose ele or pos.
          float rxy_gen_e = sqrt(pow(daughter.vx(), 2) + pow(daughter.vy(), 2));
          fRegistry.fill(HIST("Generated/hPhotonRZ"), daughter.vz(), rxy_gen_e);
          fRegistry.fill(HIST("Generated/hPhotonRxy"), daughter.vx(), daughter.vy());
          fRegistry.fill(HIST("Generated/hPhotonPhivsRxy"), daughter.phi(), rxy_gen_e);
        }
      } // end of mctrack loop per collision
    } // end of collision loop
  }

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(PCMQCMC, processQCMC, "run PCM QC in MC", false);
  PROCESS_SWITCH(PCMQCMC, processGen, "run generated information", false);
  PROCESS_SWITCH(PCMQCMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQCMC>(cfgc, TaskName{"pcm-qc-mc"})};
}
