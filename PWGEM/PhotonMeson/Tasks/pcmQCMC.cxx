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

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMEventCut.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::mcutil;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::BinnedGenPts>;
using MyMCCollision = MyMCCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

struct PCMQCMC {

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<float> maxRgen{"maxRgen", 90.f, "maximum radius for generated particles"};
  Configurable<float> margin_z_mc{"margin_z_mc", 7.0, "margin for z cut in cm for MC"};

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

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<bool> cfg_require_v0_on_wwire_ib{"cfg_require_v0_on_wwire_ib", false, "flag to select V0s on W wires ITSib"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.9, "max eta for v0 photons at PV"};
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
  } pcmcuts;

  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

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
      ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.1 GeV/c, every 0.05 GeV/c
    }
    for (int i = 2; i < 52; i++) {
      ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 5 GeV/c, every 0.1 GeV/c
    }
    for (int i = 52; i < 62; i++) {
      ptbins.emplace_back(0.5 * (i - 52) + 5.0); // from 5 to 10 GeV/c, evety 0.5 GeV/c
    }
    for (int i = 62; i < 73; i++) {
      ptbins.emplace_back(1.0 * (i - 62) + 10.0); // from 10 to 20 GeV/c, evety 1 GeV/c
    }
    const AxisSpec axis_pt{ptbins, "p_{T,#gamma} (GeV/c)"};
    const AxisSpec axis_rapidity{{0.0, +0.8, +0.9}, "rapidity |y_{#gamma}|"};
    fRegistry.add("Generated/hPt", "pT;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);
    fRegistry.add("Generated/hPtY", "Generated info", kTH2F, {axis_pt, axis_rapidity}, true);
    fRegistry.add("Generated/hPt_ConvertedPhoton", "converted photon pT;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);
    fRegistry.add("Generated/hY_ConvertedPhoton", "converted photon y;rapidity y", kTH1F, {{40, -2.0f, 2.0f}}, true);
    fRegistry.add("Generated/hPhi_ConvertedPhoton", "converted photon #varphi;#varphi (rad.)", kTH1F, {{180, 0, 2 * M_PI}}, true);
    fRegistry.add("Generated/hPhotonRxy", "conversion point in XY MC;V_{x} (cm);V_{y} (cm)", kTH2F, {{800, -100.0f, 100.0f}, {800, -100.0f, 100.0f}}, true);
    fRegistry.add("Generated/hPhotonRZ", "conversion point in RZ MC;V_{z} (cm);R_{xy} (cm)", kTH2F, {{400, -100.0f, 100.0f}, {400, 0.f, 100.0f}}, true);
    fRegistry.add("Generated/hPhotonPhivsRxy", "conversion point of #varphi vs. R_{xy} MC;#varphi (rad.);R_{xy} (cm);N_{e}", kTH2F, {{360, 0.0f, 2 * M_PI}, {400, 0, 100}}, true);

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
    fRegistry.add("V0/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
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

    fRegistry.add("V0/hPt_Photon_Candidate", "pT of photon candidate;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);                                                        // for denominator of purity
    fRegistry.add("V0/hEtaPhi_Photon_Candidate", "#eta vs. #varphi of photon candidate ;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, true); // for denominator of purity

    fRegistry.add("V0/hPt_Photon_Primary", "pT;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);                                                       // for MC efficiency
    fRegistry.add("V0/hEtaPhi_Photon_Primary", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, true); // for MC efficiency
    fRegistry.add("V0/hCosPA_Rxy_Photon_Primary", "cos PA vs. R_{xy};R_{xy} (cm);cosine pointing angle", kTH2F, {{200, 0.f, 100.f}, {100, 0.9f, 1.0f}}, true);
    fRegistry.add("V0/hXY_Photon_Primary", "X vs. Y of photon rec.;X (cm);Y (cm)", kTH2F, {{400, -100.0f, +100}, {400, -100, +100}}, true);
    fRegistry.add("V0/hXY_Photon_Primary_MC", "X vs. Y of photon gen.;X (cm);Y (cm)", kTH2F, {{400, -100.0f, +100}, {400, -100, +100}}, true);
    fRegistry.add("V0/hRZ_Photon_Primary", "R vs. Z of photon rec.;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100.0f, +100}, {200, 0, 100}}, true);
    fRegistry.add("V0/hRZ_Photon_Primary_MC", "R vs. Z of photon gen;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100.0f, +100}, {200, 0, 100}}, true);

    fRegistry.add("V0/hPt_Photon_FromWD", "pT;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);                                                       // for MC feed down correction
    fRegistry.add("V0/hEtaPhi_Photon_FromWD", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, true); // for MC feed down correction
    fRegistry.add("V0/hCosPA_Rxy_Photon_FromWD", "cos PA vs. R_{xy};R_{xy} (cm);cosine pointing angle", kTH2F, {{200, 0.f, 100.f}, {100, 0.9f, 1.0f}}, true);

    fRegistry.add("V0/hPt_Photon_hs", "pT of photon from hadronic shower in materials;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);
    fRegistry.add("V0/hEtaPhi_Photon_hs", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, true);
    fRegistry.add("V0/hCosPA_Rxy_Photon_hs", "cos PA vs. R_{xy};R_{xy} (cm);cosine pointing angle", kTH2F, {{200, 0.f, 100.f}, {100, 0.9f, 1.0f}}, true);
    fRegistry.add("V0/hXY_Photon_hs_MC", "X vs. Y of photon from hadronic shower in materials gen.;Z (cm);R_{xy} (cm)", kTH2F, {{400, -100.0f, +100}, {400, -100, 100}}, true); // production vertex of photons
    fRegistry.add("V0/hRZ_Photon_hs_MC", "R vs. Z of photon from hadronic shower in materials gen.;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100.0f, +100}, {200, 0, 100}}, true);    // production vertex of photons

    fRegistry.add("V0/hConvPoint_diffX", "conversion point diff X MC;X_{MC} (cm);X_{rec} - X_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
    fRegistry.add("V0/hConvPoint_diffY", "conversion point diff Y MC;Y_{MC} (cm);Y_{rec} - Y_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);
    fRegistry.add("V0/hConvPoint_diffZ", "conversion point diff Z MC;Z_{MC} (cm);Z_{rec} - Z_{MC} (cm)", kTH2F, {{200, -100, +100}, {100, -50.0f, 50.0f}}, true);

    fRegistry.add("V0/hPtGen_DeltaPtOverPtGen", "photon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/hPtGen_DeltaEta", "photon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/hPtGen_DeltaPhi", "photon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/hEtaRec_DeltaPtOverPtGen", "photon p_{T} resolution;#eta^{rec} of conversion point;(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{400, -2, +2}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/hEtaRec_DeltaEta", "photon #eta resolution;#eta^{rec} of conversion point;#eta^{rec} - #eta^{gen}", kTH2F, {{400, -2, +2}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0/hEtaRec_DeltaPhi", "photon #varphi resolution;#eta^{rec} of conversion point;#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{400, -2, +2}, {400, -1.0f, 1.0f}}, true);

    // v0leg info
    fRegistry.add("V0Leg/hPt", "pT;p_{T,e} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("V0Leg/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{1000, -50, 50}}, false);
    fRegistry.add("V0Leg/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
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
    fRegistry.add("V0Leg/hMeanClusterSizeITS", "mean cluster size ITS;<cluster size> on ITS #times cos(#lambda)", kTH1F, {{32, 0, 16}}, false);
    fRegistry.add("V0Leg/hXY", "X vs. Y;X (cm);Y (cm)", kTH2F, {{100, 0, 100}, {100, -50, 50}}, false);
    fRegistry.add("V0Leg/hZX", "Z vs. X;Z (cm);X (cm)", kTH2F, {{200, -100, 100}, {100, 0, 100}}, false);
    fRegistry.add("V0Leg/hZY", "Z vs. Y;Z (cm);Y (cm)", kTH2F, {{200, -100, 100}, {100, -50, 50}}, false);

    fRegistry.add("V0Leg/hPtGen_DeltaPtOverPtGen", "electron p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0Leg/hPtGen_DeltaEta", "electron #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0Leg/hPtGen_DeltaPhi", "electron #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0Leg/hEtaRec_DeltaPtOverPtGen", "electron p_{T} resolution;#eta^{rec} of conversion point;(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{400, -2, +2}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0Leg/hEtaRec_DeltaEta", "electron #eta resolution;#eta^{rec} of conversion point;#eta^{rec} - #eta^{gen}", kTH2F, {{400, -2, +2}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("V0Leg/hEtaRec_DeltaPhi", "electron #varphi resolution;#eta^{rec} of conversion point;#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{400, -2, +2}, {400, -1.0f, 1.0f}}, true);
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
  void fillV0Info(TV0 const& v0, const float weight = 1.f)
  {
    fRegistry.fill(HIST("V0/hPt"), v0.pt(), weight);
    fRegistry.fill(HIST("V0/hEtaPhi"), v0.phi(), v0.eta(), weight);
    fRegistry.fill(HIST("V0/hRadius"), v0.vz(), v0.v0radius(), weight);
    fRegistry.fill(HIST("V0/hCosPA"), v0.cospa(), weight);
    fRegistry.fill(HIST("V0/hCosPA_Rxy"), v0.v0radius(), v0.cospa(), weight);
    fRegistry.fill(HIST("V0/hPCA"), v0.pca(), weight);
    fRegistry.fill(HIST("V0/hPCA_CosPA"), v0.cospa(), v0.pca(), weight);
    fRegistry.fill(HIST("V0/hPCA_Rxy"), v0.v0radius(), v0.pca(), weight);
    fRegistry.fill(HIST("V0/hDCAxyz"), v0.dcaXYtopv(), v0.dcaZtopv(), weight);
    fRegistry.fill(HIST("V0/hAPplot"), v0.alpha(), v0.qtarm(), weight);
    fRegistry.fill(HIST("V0/hMassGamma"), v0.v0radius(), v0.mGamma(), weight);
    fRegistry.fill(HIST("V0/hGammaRxy"), v0.vx(), v0.vy(), weight);
    fRegistry.fill(HIST("V0/hKFChi2vsM"), v0.mGamma(), v0.chiSquareNDF(), weight);
    fRegistry.fill(HIST("V0/hKFChi2vsR"), v0.v0radius(), v0.chiSquareNDF(), weight);
    fRegistry.fill(HIST("V0/hKFChi2vsX"), v0.vx(), v0.chiSquareNDF(), weight);
    fRegistry.fill(HIST("V0/hKFChi2vsY"), v0.vy(), v0.chiSquareNDF(), weight);
    fRegistry.fill(HIST("V0/hKFChi2vsZ"), v0.vz(), v0.chiSquareNDF(), weight);
  }

  template <typename TLeg>
  void fillV0LegInfo(TLeg const& leg, const float weight = 1.f)
  {
    fRegistry.fill(HIST("V0Leg/hPt"), leg.pt(), weight);
    fRegistry.fill(HIST("V0Leg/hQoverPt"), leg.sign() / leg.pt(), weight);
    fRegistry.fill(HIST("V0Leg/hEtaPhi"), leg.phi(), leg.eta(), weight);
    fRegistry.fill(HIST("V0Leg/hDCAxyz"), leg.dcaXY(), leg.dcaZ(), weight);
    fRegistry.fill(HIST("V0Leg/hNclsITS"), leg.itsNCls(), weight);
    fRegistry.fill(HIST("V0Leg/hNclsTPC"), leg.tpcNClsFound(), weight);
    fRegistry.fill(HIST("V0Leg/hNcrTPC"), leg.tpcNClsCrossedRows(), weight);
    fRegistry.fill(HIST("V0Leg/hTPCNcr2Nf"), leg.tpcCrossedRowsOverFindableCls(), weight);
    fRegistry.fill(HIST("V0Leg/hTPCNcls2Nf"), leg.tpcFoundOverFindableCls(), weight);
    fRegistry.fill(HIST("V0Leg/hChi2TPC"), leg.tpcChi2NCl(), weight);
    fRegistry.fill(HIST("V0Leg/hChi2ITS"), leg.itsChi2NCl(), weight);
    fRegistry.fill(HIST("V0Leg/hITSClusterMap"), leg.itsClusterMap(), weight);
    fRegistry.fill(HIST("V0Leg/hMeanClusterSizeITS"), leg.meanClusterSizeITS() * std::cos(std::atan(leg.tgl())), weight);
    fRegistry.fill(HIST("V0Leg/hTPCdEdx"), leg.tpcInnerParam(), leg.tpcSignal(), weight);
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl(), weight);
    fRegistry.fill(HIST("V0Leg/hTPCNsigmaPi"), leg.tpcInnerParam(), leg.tpcNSigmaPi(), weight);
    fRegistry.fill(HIST("V0Leg/hXY"), leg.x(), leg.y(), weight);
    fRegistry.fill(HIST("V0Leg/hZX"), leg.z(), leg.x(), weight);
    fRegistry.fill(HIST("V0Leg/hZY"), leg.z(), leg.y(), weight);
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  void processQCMC(FilteredMyCollisions const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const&, aod::EMMCParticles const& mcparticles, MyMCCollisions const&)
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
      int nv0 = 0;
      for (auto& v0 : V0Photons_coll) {
        auto pos = v0.posTrack_as<MyMCV0Legs>();
        auto ele = v0.negTrack_as<MyMCV0Legs>();
        auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
        auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

        // LOGF(info, "posmc.isPhysicalPrimary() = %d, posmc.producedByGenerator() = %d, elemc.isPhysicalPrimary() = %d, elemc.producedByGenerator() = %d", posmc.isPhysicalPrimary(), posmc.producedByGenerator(), elemc.isPhysicalPrimary(), elemc.producedByGenerator());

        if (fV0PhotonCut.IsSelected<MyMCV0Legs>(v0)) {
          fRegistry.fill(HIST("V0/hPt_Photon_Candidate"), v0.pt());
          fRegistry.fill(HIST("V0/hEtaPhi_Photon_Candidate"), v0.phi(), v0.eta());

          int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
          if (photonid < 0) {
            continue;
          }
          auto mcphoton = mcparticles.iteratorAt(photonid);
          nv0++; // number of true photons regardless of primary or feeddown or hadronic shower

          if (IsConversionPointInAcceptance(mcphoton, maxRgen, pcmcuts.cfg_max_eta_v0, margin_z_mc, mcparticles)) {
            if (mcphoton.isPhysicalPrimary() || mcphoton.producedByGenerator()) {
              float rxy_rec = sqrt(v0.vx() * v0.vx() + v0.vy() * v0.vy());
              float rxy_mc = sqrt(posmc.vx() * posmc.vx() + posmc.vy() * posmc.vy());
              float eta_cp = std::atanh(v0.vz() / sqrt(pow(v0.vx(), 2) + pow(v0.vy(), 2) + pow(v0.vz(), 2)));
              fillV0Info(v0, 1.f);
              fRegistry.fill(HIST("V0/hPt_Photon_Primary"), v0.pt());
              fRegistry.fill(HIST("V0/hEtaPhi_Photon_Primary"), v0.phi(), v0.eta());
              fRegistry.fill(HIST("V0/hXY_Photon_Primary"), v0.vx(), v0.vy());
              fRegistry.fill(HIST("V0/hXY_Photon_Primary_MC"), posmc.vx(), posmc.vy());
              fRegistry.fill(HIST("V0/hRZ_Photon_Primary"), v0.vz(), rxy_rec);
              fRegistry.fill(HIST("V0/hRZ_Photon_Primary_MC"), posmc.vz(), rxy_mc);
              fRegistry.fill(HIST("V0/hCosPA_Rxy_Photon_Primary"), v0.v0radius(), v0.cospa());
              fRegistry.fill(HIST("V0/hPtGen_DeltaPtOverPtGen"), mcphoton.pt(), (v0.pt() - mcphoton.pt()) / mcphoton.pt());
              fRegistry.fill(HIST("V0/hPtGen_DeltaEta"), mcphoton.pt(), v0.eta() - mcphoton.eta());
              fRegistry.fill(HIST("V0/hPtGen_DeltaPhi"), mcphoton.pt(), v0.phi() - mcphoton.phi());
              fRegistry.fill(HIST("V0/hEtaRec_DeltaPtOverPtGen"), eta_cp, (v0.pt() - mcphoton.pt()) / mcphoton.pt());
              fRegistry.fill(HIST("V0/hEtaRec_DeltaEta"), eta_cp, v0.eta() - mcphoton.eta());
              fRegistry.fill(HIST("V0/hEtaRec_DeltaPhi"), eta_cp, v0.phi() - mcphoton.phi());
              fRegistry.fill(HIST("V0/hConvPoint_diffX"), elemc.vx(), v0.vx() - elemc.vx());
              fRegistry.fill(HIST("V0/hConvPoint_diffY"), elemc.vy(), v0.vy() - elemc.vy());
              fRegistry.fill(HIST("V0/hConvPoint_diffZ"), elemc.vz(), v0.vz() - elemc.vz());

              for (auto& leg : {pos, ele}) {
                fillV0LegInfo(leg, 1.f);
                auto mcleg = leg.template emmcparticle_as<aod::EMMCParticles>();
                fRegistry.fill(HIST("V0Leg/hPtGen_DeltaPtOverPtGen"), mcleg.pt(), (leg.pt() - mcleg.pt()) / mcleg.pt());
                fRegistry.fill(HIST("V0Leg/hPtGen_DeltaEta"), mcleg.pt(), leg.eta() - mcleg.eta());
                fRegistry.fill(HIST("V0Leg/hPtGen_DeltaPhi"), mcleg.pt(), leg.phi() - mcleg.phi());
                fRegistry.fill(HIST("V0Leg/hEtaRec_DeltaPtOverPtGen"), eta_cp, (leg.pt() - mcleg.pt()) / mcleg.pt());
                fRegistry.fill(HIST("V0Leg/hEtaRec_DeltaEta"), eta_cp, leg.eta() - mcleg.eta());
                fRegistry.fill(HIST("V0Leg/hEtaRec_DeltaPhi"), eta_cp, leg.phi() - mcleg.phi());
              }
            } else if (IsFromWD(mcphoton.template emmcevent_as<MyMCCollisions>(), mcphoton, mcparticles) > 0) {
              fRegistry.fill(HIST("V0/hPt_Photon_FromWD"), v0.pt());
              fRegistry.fill(HIST("V0/hEtaPhi_Photon_FromWD"), v0.phi(), v0.eta());
              fRegistry.fill(HIST("V0/hCosPA_Rxy_Photon_FromWD"), v0.v0radius(), v0.cospa());
            } else {
              // LOGF(info, "mcphoton.vx() = %f, mcphoton.vy() = %f, mcphoton.vz() = %f, mother_pdg = %d", mcphoton.vx(), mcphoton.vy(), mcphoton.vz(), mother_pdg);
              float rxy_photon_hs = sqrt(mcphoton.vx() * mcphoton.vx() + mcphoton.vy() * mcphoton.vy());
              fRegistry.fill(HIST("V0/hPt_Photon_hs"), v0.pt());
              fRegistry.fill(HIST("V0/hEtaPhi_Photon_hs"), v0.phi(), v0.eta());
              fRegistry.fill(HIST("V0/hXY_Photon_hs_MC"), mcphoton.vx(), mcphoton.vy());
              fRegistry.fill(HIST("V0/hRZ_Photon_hs_MC"), mcphoton.vz(), rxy_photon_hs);
              fRegistry.fill(HIST("V0/hCosPA_Rxy_Photon_hs"), v0.v0radius(), v0.cospa());
            }
          }
        }
      } // end of v0 loop
      fRegistry.fill(HIST("V0/hNgamma"), nv0);
    } // end of collision loop
  }   // end of process

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

        if (abs(mctrack.pdgCode()) == 22 && (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) && IsConversionPointInAcceptance(mctrack, maxRgen, pcmcuts.cfg_max_eta_v0, margin_z_mc, mcparticles)) {
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
    }   // end of collision loop
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
