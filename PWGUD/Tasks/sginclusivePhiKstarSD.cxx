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
/// \file sginclusivePhiKstarSD.cxx
/// \brief Single Gap Event Analyzer for phi and Kstar
/// \author Sandeep Dudi <sandeep.dudi3@gmail.com>, Subhadeep Mandal <subhadeep.mandal@cern.ch>
/// \since  May 2024

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TPDGCode.h"
#include <TMath.h>
#include <TString.h>

#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

// GapSide enum
using o2::aod::sgselector::DoubleGap;
using o2::aod::sgselector::SingleGapA;
using o2::aod::sgselector::SingleGapC;

struct SginclusivePhiKstarSD {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rQA{"QA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 100., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> fddaCut{"fddaCut", 10000., "FDDA threshold"};
  Configurable<float> fddcCut{"fddcCut", 10000., "FDDC threshold"};
  Configurable<float> zdcCut{"zdcCut", 0., "ZDC threshold"};
  Configurable<float> vzCut{"vzCut", 10., "Vz position"};
  Configurable<bool> useOccCut{"useOccCut", false, "Turn on/off Occupancy cut"};
  Configurable<float> confgOccCut{"confgOccCut", 1000., "Occupancy cut"};
  Configurable<bool> useHadronicRateCut{"useHadronicRateCut", false, "Turn on/off hadronicRate cut"};
  Configurable<float> confgHadronicRateMax{"confgHadronicRateMax", 1000., "Maximum hadronicRate cut"};
  Configurable<float> confgHadronicRateMin{"confgHadronicRateMin", 0., "Minimum hadronicRate cut"};
  Configurable<bool> useTrs{"useTrs", false, "kNoCollInTimeRangeStandard cut"};
  Configurable<bool> useTrofs{"useTrofs", false, "kNoCollInRofStandard cut"};
  Configurable<bool> useHmpr{"useHmpr", false, "kNoHighMultCollInPrevRof cut"};
  Configurable<bool> useTfb{"useTfb", false, "kNoTimeFrameBorder cut"};
  Configurable<bool> useItsrofb{"useItsrofb", false, "kNoITSROFrameBorder cut"};
  Configurable<bool> useSbp{"useSbp", false, "kNoSameBunchPileup cut"};
  Configurable<bool> useZvtxftovpv{"useZvtxftovpv", false, "kIsGoodZvtxFT0vsPV cut"};
  Configurable<bool> useVtxItsTpc{"useVtxItsTpc", false, "kIsVertexITSTPC cut"};
  Configurable<bool> usenumContrib{"usenumContrib", false, "numContrib cut for event selection"};
  Configurable<int> upcflag{"upcflag", -1, "upc run selection, -1 = off, 0 = std, 1 = upc"};
  Configurable<int> cutRCTflag{"cutRCTflag", 0, {"0 = off, 1 = CBT, 2 = CBT+ZDC, 3 = CBThadron, 4 = CBThadron+ZDC"}};

  // Track Selections
  Configurable<float> pvCut{"pvCut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcazCut{"dcazCut", 2.0, "dcaZ cut"};
  Configurable<float> dcaxyCut{"dcaxyCut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "Max itsChi2NCl"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> ptCut{"ptCut", 0.15, "Track pt cut"};
  Configurable<float> pt1{"pt1", 0.3, "pid selection pt1"};
  Configurable<float> pt2{"pt2", 0.4, "pid selection pt2"};
  Configurable<float> pt3{"pt3", 0.5, "pid selection pt3"};

  Configurable<bool> useMultCut{"useMultCut", false, "Multipicity cut on good tracks"};

  Configurable<bool> rapiditycut{"rapiditycut", true, "Rapidity Cut"};
  Configurable<float> rapiditycutvalue{"rapiditycutvalue", 0.5, "Rapidity Cut value"};

  Configurable<float> nsigmaTpcCut1{"nsigmaTpcCut1", 3.0, "nsigma tpc cut1"};
  Configurable<float> nsigmaTpcCut2{"nsigmaTpcCut2", 3.0, "nsigma tpc cut2"};
  Configurable<float> nsigmaTpcCut3{"nsigmaTpcCut3", 3.0, "nsigma tpc cut3"};
  Configurable<float> nsigmaTpcCut{"nsigmaTpcCut", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmaTofCut{"nsigmaTofCut", 9.0, "nsigma tpc+tof cut"};
  Configurable<float> nsigmaTofCut1{"nsigmaTofCut1", 3.0, "nsigma tof cut"};
  Configurable<float> pionNsigmaCut{"pionNsigmaCut", 3.0, "nsigma tpc cut for kaon"};

  Configurable<int> mintrack{"mintrack", 1, "min track"};
  Configurable<int> maxtrack{"maxtrack", 150, "max track"};
  Configurable<bool> useTof{"useTof", true, "TOF PID"};
  Configurable<bool> ccut{"ccut", true, "TPC + TOF PID"};
  Configurable<bool> kaoncut{"kaoncut", false, " kaon slection cut for kstar "};

  Configurable<bool> qa{"qa", false, "QA plots for Data (turn qaMC to 0)"};
  Configurable<bool> qaMC{"qaMC", false, "QA plots for MC (turn qa for data to 0)"};
  Configurable<bool> exclusive{"exclusive", false, "for double gap side "};

  Configurable<bool> phi{"phi", true, ""};
  Configurable<bool> rho{"rho", true, ""};
  Configurable<bool> kstar{"kstar", true, ""};
  Configurable<bool> fourpion{"fourpion", true, ""};
  Configurable<bool> mc{"mc", true, ""};
  Configurable<int> gapsideMC{"gapsideMC", 1, "gapside MC"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * o2::constants::math::PI / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * o2::constants::math::PI / 6.0, "Maximum of rotation"};
  //
  Configurable<bool> reconstruction{"reconstruction", true, ""};
  Configurable<int> generatedId{"generatedId", 40, "40 = PhiA, 44 PhiC, 41 = K*0A, 45 = K*0C"};

  ConfigurableAxis axisphimass{"axisphimass", {220, 0.98, 1.2}, ""};
  ConfigurableAxis axiskstarmass{"axiskstarmass", {400, 0.0, 2.0}, ""};
  ConfigurableAxis axisrhomass{"axisrhomass", {200, 1.0, 2.0}, ""};
  ConfigurableAxis axispt{"axispt", {200, 0.0, 20.0}, ""};
  ConfigurableAxis axisrapdity{"axisrapdity", {40, -2.0, 2.0}, ""};

  int numTwoTracks = 2;
  int numFourTracks = 4;

  void init(InitContext const& context)
  {
    registry.add("hEventCutFlow", "No. of events after event cuts", kTH1F, {{20, 0, 20}});
    std::shared_ptr<TH1> hCutFlow = registry.get<TH1>(HIST("hEventCutFlow"));

    auto check = [](bool enabled) { return enabled ? "" : " #otimes"; }; // check if a cut is enabled and put #otimes if not enabled beside that label

    std::vector<std::string> eveCutLabels = {
      "All Events",
      "Gapside (0 to 2)",
      Form("|Vz| < %.1f", vzCut.value),
      Form("Occupancy < %.0f%s", confgOccCut.value, check(useOccCut.value)),
      Form("%.1e < Hadronic Rate < %.1e%s", confgHadronicRateMin.value, confgHadronicRateMax.value, check(useHadronicRateCut.value)),
      std::string("kNoCollInTimeRangeStandard") + check(useTrs.value),
      std::string("kNoCollInRofStandard") + check(useTrofs.value),
      std::string("kNoHighMultCollInPrevRof") + check(useHmpr.value),
      std::string("kNoTimeFrameBorder") + check(useTfb.value),
      std::string("kNoITSROFrameBorder") + check(useItsrofb.value),
      std::string("kNoSameBunchPileup") + check(useSbp.value),
      std::string("kIsGoodZvtxFT0vsPV") + check(useZvtxftovpv.value),
      std::string("kIsVertexITSTPC") + check(useVtxItsTpc.value),
      Form("RCTFlag = %d%s", cutRCTflag.value, check(cutRCTflag.value > 0)),
      Form("upcFlag = %d%s", upcflag.value, check(upcflag.value != -1)),
      Form("%d < numContrib < %d%s", mintrack.value, maxtrack.value, check(usenumContrib.value))};

    for (size_t i = 0; i < eveCutLabels.size(); ++i) {
      hCutFlow->GetXaxis()->SetBinLabel(i + 1, eveCutLabels[i].c_str());
    }

    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("nPVContributors_data", "Multiplicity_dist_before track cut gap A", kTH1F, {{110, 0, 110}});
    registry.add("nPVContributors_data_1", "Multiplicity_dist_before track cut gap C", kTH1F, {{110, 0, 110}});

    registry.add("hRotation", "hRotation", kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});

    if (phi) {
      registry.add("os_KK_pT_0", "pt kaon pair", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_pT_1", "pt kaon pair", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_pT_2", "pt kaon pair", kTH3F, {axisphimass, axisrapdity, axispt});

      registry.add("os_KK_lsMM_pT_0", "kaon pair Negative like sign", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_lsPP_pT_0", "kaon pair Positive like sign", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_lsMM_pT_1", "kaon pair Negative like sign", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_lsPP_pT_1", "kaon pair Positive like sign", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_ls_pT_2", "kaon pair like sign", kTH3F, {axisphimass, axisrapdity, axispt});

      registry.add("os_KK_mix_pT_0", "kaon pair mix event", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_mix_pT_1", "kaon pair mix event", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_mix_pT_2", "kaon pair mix event", kTH3F, {axisphimass, axisrapdity, axispt});

      registry.add("os_KK_rot_pT_0", "kaon pair rotional event", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_rot_pT_1", "kaon pair rotional event", kTH3F, {axisphimass, axisrapdity, axispt});
      registry.add("os_KK_rot_pT_2", "kaon pair rotional event", kTH3F, {axisphimass, axisrapdity, axispt});
    }
    if (rho) {
      registry.add("os_pp_pT_0", "pt pion pair", kTH3F, {axisrhomass, axisrapdity, axispt});
      registry.add("os_pp_pT_1", "pt pion pair", kTH3F, {axisrhomass, axisrapdity, axispt});
      registry.add("os_pp_pT_2", "pt pion pair", kTH3F, {axisrhomass, axisrapdity, axispt});
      registry.add("os_pp_ls_pT_0", "pion pair like sign", kTH3F, {axisrhomass, axisrapdity, axispt});
      registry.add("os_pp_ls_pT_1", "pion pair like sign", kTH3F, {axisrhomass, axisrapdity, axispt});
      registry.add("os_pp_ls_pT_2", "pion pair like sign", kTH3F, {axisrhomass, axisrapdity, axispt});
    }
    if (kstar) {
      registry.add("os_pk_pT_0", "pion-kaon pair", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_pT_1", "pion-kaon pair", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_pT_2", "pion-kaon pair", kTH3F, {axiskstarmass, axisrapdity, axispt});

      registry.add("os_pk_mix_pT_0", "pion-kaon mix pair", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_mix_pT_1", "pion-kaon mix pair", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_mix_pT_2", "pion-kaon mix pair", kTH3F, {axiskstarmass, axisrapdity, axispt});

      registry.add("os_pk_rot_pT_0", "pion-kaon rotional pair", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_rot_pT_1", "pion-kaon rotional pair", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_rot_pT_2", "pion-kaon rotional pair", kTH3F, {axiskstarmass, axisrapdity, axispt});

      registry.add("os_pk_lsMM_pT_0", "pion-kaon pair Negative like sign", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_lsPP_pT_0", "pion-kaon pair Positive like sign", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_lsMM_pT_1", "pion-kaon pair Negative like sign", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_lsPP_pT_1", "pion-kaon pair Positive like sign", kTH3F, {axiskstarmass, axisrapdity, axispt});
      registry.add("os_pk_ls_pT_2", "pion-kaon like sign", kTH3F, {axiskstarmass, axisrapdity, axispt});
    }
    // qa plots
    if (qa) {
      // Occupancy
      rQA.add("hOcc_before", "Occupancy distribution before event cuts", kTH1F, {{1000, 0, 50000}});
      rQA.add("hOcc_after", "Occupancy distribution after all event cuts", kTH1F, {{1000, 0, 10000}});

      // DCA
      rQA.add("hDcaxy_all_before", "DCAxy Distribution of all tracks before track selection; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
      rQA.add("hDcaz_all_before", "DCAz Distribution of all tracks before track selection; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

      rQA.add("hDcaxy_all_after", "DCAxy Distribution of all tracks after track selection; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
      rQA.add("hDcaz_all_after", "DCAz Distribution of all tracks after track selection; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

      rQA.add("hDcaxy_pi", "DCAxy Distribution of selected pions; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
      rQA.add("hDcaz_pi", "DCAz Distribution of selected pions; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

      rQA.add("hDcaxy_ka", "DCAxy Distribution of selected kaons; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
      rQA.add("hDcaz_ka", "DCAz Distribution of selected kaons; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

      // Vx, Vy, Vz
      rQA.add("hVertexX", "Vertex X distribution; Vertex X (cm); Counts", kTH1F, {{400, -0.1, 0.1}});
      rQA.add("hVertexY", "Vertex Y distribution; Vertex Y (cm); Counts", kTH1F, {{200, -0.05, 0.05}});
      rQA.add("hVertexZ", "VertexZ distribution; Vertex Z (cm); Counts", kTH1F, {{600, -15.0, 15.0}});

      // TPC, TOF PID
      rQA.add("tpc_dedx", "p vs dE/dx of all particles; #it{p} (GeV/#it{c}); TPC dE/dx (a.u.)", kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      rQA.add("tpc_dedx_kaon", "p vs dE/dx of selected kaons; #it{p} (GeV/#it{c}); TPC dE/dx (a.u.)", kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      rQA.add("tpc_dedx_pion", "p vs dE/dx of selected pions; #it{p} (GeV/#it{c}); TPC dE/dx (a.u.)", kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      rQA.add("tof_beta", "p vs #beta of all particles; #it{p} (GeV/#it{c}); TOF #beta", kTH2F, {{500, 0.0, 10.0}, {500, 0.0, 1.0}});
      rQA.add("tof_beta_kaon", "p vs #beta of selected kaons; #it{p} (GeV/#it{c}); TOF #beta", kTH2F, {{500, 0.0, 10.0}, {500, 0.0, 1.0}});
      rQA.add("tof_beta_pion", "p vs #beta of selected pions; #it{p} (GeV/#it{c}); TOF #beta", kTH2F, {{500, 0.0, 10.0}, {500, 0.0, 1.0}});

      rQA.add("tpc_nsigma_kaon_all", "Kaon n#sigma_{TPC} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      rQA.add("tpc_nsigma_pion_all", "Pion n#sigma_{TPC} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      rQA.add("tpc_nsigma_kaon", "Kaon n#sigma_{TPC} of selected kaons; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      rQA.add("tpc_nsigma_pion", "Pion n#sigma_{TPC} of selected pions; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

      rQA.add("tof_nsigma_kaon_all", "Kaon n#sigma_{TOF} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      rQA.add("tof_nsigma_pion_all", "Pion n#sigma_{TOF} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      rQA.add("tof_nsigma_kaon", "Kaon n#sigma_{TOF} of selected kaons; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      rQA.add("tof_nsigma_pion", "Pion n#sigma_{TPC} of selected pions; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

      rQA.add("tpc_tof_nsigma_kaon", "n#sigma TPC vs TOF; n#sigma_{TPC}^{K}; n#sigma_{TOF}^{K}", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});
      rQA.add("tpc_tof_nsigma_pion", "n#sigma TPC vs TOF; n#sigma_{TPC}^{#pi}; n#sigma_{TOF}^{#pi}", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});

      // Rapidity, pseudorapisdity
      rQA.add("hEta_all_after", "Pseudorapidity of all tracks after track selection; #eta; Counts", kTH1F, {{400, -1.0, 1.0}});

      rQA.add("hEta_ka", "Pseudorapidity of selected Kaons; #eta; Counts", kTH1F, {{400, -1.0, 1.0}});
      rQA.add("hRap_ka", "Rapidity of selected Kaons; y; Counts", kTH1F, {{400, -1.0, 1.0}});

      rQA.add("hEta_pi", "Pseudorapidity of selected Pions; #eta; Counts", kTH1F, {{400, -1.0, 1.0}});
      rQA.add("hRap_pi", "Rapidity of selected Pions; y; Counts", kTH1F, {{400, -1.0, 1.0}});

      // Detector Signals
      rQA.add("FT0A_2", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
      rQA.add("FT0A_0", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
      rQA.add("FT0A_1", "T0A amplitude", kTH1F, {{20000, 0.0, 20000.0}});
      rQA.add("FT0C_2", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
      rQA.add("FT0C_0", "T0C amplitude", kTH1F, {{20000, 0.0, 20000.0}});
      rQA.add("FT0C_1", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
      rQA.add("ZDC_A_2", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      rQA.add("ZDC_A_0", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      rQA.add("ZDC_A_1", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      rQA.add("ZDC_C_2", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      rQA.add("ZDC_C_0", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      rQA.add("ZDC_C_1", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      rQA.add("V0A_2", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
      rQA.add("V0A_0", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
      rQA.add("V0A_1", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
    }
    registry.add("gap_mult0", "Mult 0", kTH1F, {{100, 0.0, 100.0}});
    registry.add("gap_mult1", "Mult 1", kTH1F, {{100, 0.0, 100.0}});
    registry.add("gap_mult2", "Mult 2", kTH1F, {{100, 0.0, 100.0}});

    registry.add("mult_0", "mult0", kTH1F, {{150, 0, 150}});
    registry.add("mult_1", "mult1", kTH1F, {{150, 0, 150}});
    registry.add("mult_2", "mult2", kTH1F, {{150, 0, 150}});

    if (fourpion) {
      registry.add("os_pppp_pT_2", "4 pion pair", kTH3F, {{800, 0.5, 4.5}, {250, 0.0, 5.0}, {30, -1.5, 1.5}});
      registry.add("os_pppp_pT_2_ls", "4 pion pair", kTH3F, {{800, 0.5, 4.5}, {250, 0.0, 5.0}, {30, -1.5, 1.5}});
      registry.add("os_pp_vs_pp_mass", "pair1 vd pair2 ", kTH2F, {{800, 0.5, 4.5}, {800, 0.5, 4.5}});
      registry.add("os_pp_vs_pp_pt", "pair1 pt vs pair2 pt", kTH2F, {{250, 0.0, 5.0}, {250, 0.0, 5.0}});
      registry.add("os_pp_vs_pp_mass1", "pair3 vd pair4 ", kTH2F, {{800, 0.5, 4.5}, {800, 0.5, 4.5}});
      registry.add("os_pp_vs_pp_pt1", "pair3 pt vs pair4 pt", kTH2F, {{250, 0.0, 5.0}, {250, 0.0, 5.0}});
      registry.add("phi_dis", "phi_dis", kTH1F, {{360, 0, 6.28}});
      registry.add("costheta_dis", "costheta_dis", kTH1F, {{40, -1.0, 1.0}});
      registry.add("costheta_vs_phi", "costheta_vs_phi", kTH2F, {{40, -1.0, 1.0}, {360, 0.0, 6.28}});
      registry.add("phi_dis1", "phi_dis1", kTH1F, {{360, 0, 6.28}});
      registry.add("costheta_dis1", "costheta_dis1", kTH1F, {{40, -1.0, 1.0}});
      registry.add("costheta_vs_phi1", "costheta_vs_phi1", kTH2F, {{40, -1.0, 1.0}, {360, 0.0, 6.28}});
    }
    if (mc) {
      // add histograms for the different process functions
      if (context.mOptions.get<bool>("processMCTruth")) {
        registry.add("MC/Stat", "Count generated events; ; Entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
        registry.add("MC/recCols", "Number of reconstructed collisions; Number of reconstructed collisions; Entries", {HistType::kTH1F, {{31, -0.5, 30.5}}});
        registry.add("MC/nParts", "Number of McParticles per collision; Number of McParticles; Entries", {HistType::kTH1F, {{1001, -0.5, 1000.5}}});
        registry.add("MC/nRecTracks", "Number of reconstructed tracks per McParticle; Number of reconstructed tracks per McParticle; Entries", {HistType::kTH1F, {{11, -0.5, 10.5}}});
        registry.add("MC/genEtaPt", "Generated events; eta (1); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/genRap", "Generated events; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/genMPt", "Generated events; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{220, 0.98, 1.2}, {200, 0.0, 10.0}}});
        registry.add("MC/genM", "Generated events; Mass (GeV/#it{c}^2)", {HistType::kTH1F, {{220, 0.98, 1.2}}});
        registry.add("MC/genM_1", "Generated events; Mass (GeV/#it{c}^2)", {HistType::kTH1F, {{220, 0.98, 1.2}}});

        registry.add("MC/accMPtRap_phi_G", "Generated Phi; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_phi_G1", "Generated Phi1; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_phi_G2", "Generated Phi2; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});

        registry.add("MC/accEtaPt", "Generated events in acceptance; eta (1); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/accRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/accMPt", "Generated events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}}});
        registry.add("MC/accMPtRap", "Generated events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accM", "Generated events in acceptance; Mass (GeV/#it{c}^2)", {HistType::kTH1F, {{220, 0.98, 1.20}}});
        registry.add("MC/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/selMPt", "Selected events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
        registry.add("MC/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});
        // K*0
        registry.add("MC/accMPtRap_kstar_G", "Generated K*0; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_kstar_G1", "Generated K*0 from kp; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_kstar_G2", "Generated K*0 from kp2; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/genEtaPt_k", "Generated events; eta (1); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/genRap_k", "Generated events; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/genMPt_k", "Generated events; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("MC/genM_k", "Generated events; Mass (GeV/#it{c}^2)", {HistType::kTH1F, {{400, 0., 2.0}}});
        registry.add("MC/genM_1_k", "Generated events; Mass (GeV/#it{c}^2)", {HistType::kTH1F, {{400, 0., 2.0}}});

        registry.add("MC/accEtaPt_k", "Generated events in acceptance; eta (1); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/accRap_k", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/accMPt_k", "Generated events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("MC/accMPtRap_k", "Generated events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accM_k", "Generated events in acceptance; Mass (GeV/#it{c}^2)", {HistType::kTH1F, {{400, 0., 2.0}}});
      }
      if (context.mOptions.get<bool>("processReco")) {

        registry.add("Reco/hEventCutFlowMC", "No. of events after event cuts in MC", kTH1F, {{20, 0, 20}});
        std::shared_ptr<TH1> hCutFlowMC = registry.get<TH1>(HIST("Reco/hEventCutFlowMC"));

        std::vector<std::string> eveCutLabelsMC = {
          "All Events",
          "has_udMcCollision",
          Form("generatorsID = %d", generatedId.value),
          Form("GapsideMC = %d", gapsideMC.value),
          Form("|Vz| < %.1f", vzCut.value),
          Form("Occupancy < %.0f%s", confgOccCut.value, check(useOccCut.value)),
          Form("%.1e < Hadronic Rate < %.1e%s", confgHadronicRateMin.value, confgHadronicRateMax.value, check(useHadronicRateCut.value)),
          std::string("kNoCollInTimeRangeStandard") + check(useTrs.value),
          std::string("kNoCollInRofStandard") + check(useTrofs.value),
          std::string("kNoHighMultCollInPrevRof") + check(useHmpr.value),
          std::string("kNoTimeFrameBorder") + check(useTfb.value),
          std::string("kNoITSROFrameBorder") + check(useItsrofb.value),
          std::string("kNoSameBunchPileup") + check(useSbp.value),
          std::string("kIsGoodZvtxFT0vsPV") + check(useZvtxftovpv.value),
          std::string("kIsVertexITSTPC") + check(useVtxItsTpc.value),
          Form("RCTFlag = %d%s", cutRCTflag.value, check(cutRCTflag.value > 0)),
          Form("upcFlag = %d%s", upcflag.value, check(upcflag.value != -1)),
          Form("%d < numContrib < %d%s", mintrack.value, maxtrack.value, check(usenumContrib.value))};

        for (size_t i = 0; i < eveCutLabelsMC.size(); ++i) {
          hCutFlowMC->GetXaxis()->SetBinLabel(i + 1, eveCutLabelsMC[i].c_str());
        }

        registry.add("Reco/Stat", "Count reconstruted events; ; Entries", {HistType::kTH1F, {{5, -0.5, 4.5}}});
        registry.add("Reco/nPVContributors", "Number of PV contributors per collision; Number of PV contributors; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
        registry.add("Reco/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/selMPt", "Reconstructed events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}}});
        registry.add("Reco/selMPtRap", "Reconstructed events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selMPtRap_gen", "Reconstructed(gen) events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_phi_T", "Reconstrcted Phi; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_kstar_T", "Reconstructed K*0; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});

        registry.add("Reco/selPt", "Reconstructed events in acceptance;#it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{200, 0.0, 10.0}}});
        registry.add("Reco/selM", "Reconstructed events in acceptance; Mass (GeV/#it{c}^2); ", {HistType::kTH1F, {{220, 0.98, 1.20}}});
        registry.add("Reco/mcEtaPt", "Generated events in acceptance; eta (1); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("Reco/mcRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/mcMPt", "Generated events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
        registry.add("Reco/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});

        registry.add("Reco/selRap_k", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/selMPt_k", "Reconstructed events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("Reco/selMPtRap_k", "Reconstructed events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selMPtRap_k_gen", "Reconstructed(gen) events in acceptance; Mass (GeV/#it{c}^2); #it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selPt_k", "Reconstructed events in acceptance;#it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{200, 0.0, 10.0}}});
        registry.add("Reco/selM_k", "Reconstructed events in acceptance; Mass (GeV/#it{c}^2); ", {HistType::kTH1F, {{400, 0., 2.0}}});
        registry.add("Reco/selM_k_K", "Reconstructed events in acceptance; Mass (GeV/#it{c}^2); ", {HistType::kTH1F, {{400, 0., 2.0}}});

        registry.add("Reco/nTracks", "Number of reconstructed tracks per collision; Number of reconstructed tracks; Entries", {HistType::kTH1F, {{101, -0.5, 100.5}}});
        registry.add("Reco/treta_k", "track kaon eta", {HistType::kTH1F, {{200, -5.0, 5.0}}});
        registry.add("Reco/trpt_k", "rec kaon track pt", {HistType::kTH1F, {{200, 0.0, 10.0}}});
        registry.add("Reco/trpt", "rec track pt", {HistType::kTH1F, {{200, 0.0, 10.0}}});

        registry.add("Reco/tr_dcaz_1", "dcaz-", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_dcaxy_1", "dcaxy-", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_chi2ncl_1", "chi2ncl-", {HistType::kTH1F, {{100, 0.0, 100.0}}});
        registry.add("Reco/tr_tpcnclfind_1", "tpcnclfind-", {HistType::kTH1F, {{300, 0.0, 300.0}}});
        registry.add("Reco/tr_itsChi2NCl_1", "itsChi2NCl-", {HistType::kTH1F, {{200, 0.0, 200.0}}});
        registry.add("Reco/tr_Eta_1", " eta -", {HistType::kTH1F, {{300, 1.5, 1.5}}});

        registry.add("Reco/tr_dcaz_2", "dcaz", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_dcaxy_2", "dcaxy", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_chi2ncl_2", "chi2ncl", {HistType::kTH1F, {{100, 0.0, 100.0}}});
        registry.add("Reco/tr_tpcnclfind_2", "tpcnclfind", {HistType::kTH1F, {{300, 0.0, 300.0}}});
        registry.add("Reco/tr_itsChi2NCl_2", "itsChi2NCl", {HistType::kTH1F, {{200, 0.0, 200.0}}});

        // QA plots
        if (qaMC) {
          // Occupancy
          rQA.add("hOcc_before_mc", "Occupancy distribution before event cuts", kTH1F, {{1000, 0, 50000}});
          rQA.add("hOcc_after_mc", "Occupancy distribution after all event cuts", kTH1F, {{1000, 0, 10000}});

          // DCA
          rQA.add("hDcaxy_all_before_mc", "DCAxy Distribution of all tracks before track selection; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
          rQA.add("hDcaz_all_before_mc", "DCAz Distribution of all tracks before track selection; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

          rQA.add("hDcaxy_all_after_mc", "DCAxy Distribution of all tracks after track selection; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
          rQA.add("hDcaz_all_after_mc", "DCAz Distribution of all tracks after track selection; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

          rQA.add("hDcaxy_pi_mc", "DCAxy Distribution of selected pions; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
          rQA.add("hDcaz_pi_mc", "DCAz Distribution of selected pions; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

          rQA.add("hDcaxy_ka_mc", "DCAxy Distribution of selected kaons; DCAxy (cm); Counts", kTH1F, {{400, -0.2, 0.2}});
          rQA.add("hDcaz_ka_mc", "DCAz Distribution of selected kaons; DCAz (cm); Counts", kTH1F, {{400, -0.2, 0.2}});

          // Vx, Vy, Vz
          rQA.add("hVertexX_mc", "Vertex X distribution; Vertex X (cm); Counts", kTH1F, {{400, -0.1, 0.1}});
          rQA.add("hVertexY_mc", "Vertex Y distribution; Vertex Y (cm); Counts", kTH1F, {{200, -0.05, 0.05}});
          rQA.add("hVertexZ_mc", "VertexZ distribution; Vertex Z (cm); Counts", kTH1F, {{600, -15.0, 15.0}});

          // TPC, TOF PID
          rQA.add("tpc_dedx_mc", "p vs dE/dx of all particles; #it{p} (GeV/#it{c}); TPC dE/dx (a.u.)", kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}});
          rQA.add("tpc_dedx_kaon_mc", "p vs dE/dx of selected kaons; #it{p} (GeV/#it{c}); TPC dE/dx (a.u.)", kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}});
          rQA.add("tpc_dedx_pion_mc", "p vs dE/dx of selected pions; #it{p} (GeV/#it{c}); TPC dE/dx (a.u.)", kTH2F, {{500, 0.0, 10.0}, {5000, 0.0, 5000.0}});
          rQA.add("tof_beta_mc", "p vs #beta of all particles; #it{p} (GeV/#it{c}); TOF #beta", kTH2F, {{500, 0.0, 10.0}, {500, 0.0, 1.0}});
          rQA.add("tof_beta_kaon_mc", "p vs #beta of selected kaons; #it{p} (GeV/#it{c}); TOF #beta", kTH2F, {{500, 0.0, 10.0}, {500, 0.0, 1.0}});
          rQA.add("tof_beta_pion_mc", "p vs #beta of selected pions; #it{p} (GeV/#it{c}); TOF #beta", kTH2F, {{500, 0.0, 10.0}, {500, 0.0, 1.0}});

          rQA.add("tpc_nsigma_kaon_all_mc", "Kaon n#sigma_{TPC} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
          rQA.add("tpc_nsigma_pion_all_mc", "Pion n#sigma_{TPC} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
          rQA.add("tpc_nsigma_kaon_mc", "Kaon n#sigma_{TPC} of selected kaons; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
          rQA.add("tpc_nsigma_pion_mc", "Pion n#sigma_{TPC} of selected pions; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

          rQA.add("tof_nsigma_kaon_all_mc", "Kaon n#sigma_{TOF} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
          rQA.add("tof_nsigma_pion_all_mc", "Pion n#sigma_{TOF} of all tracks; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
          rQA.add("tof_nsigma_kaon_mc", "Kaon n#sigma_{TOF} of selected kaons; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
          rQA.add("tof_nsigma_pion_mc", "Pion n#sigma_{TPC} of selected pions; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

          rQA.add("tpc_tof_nsigma_kaon_mc", "n#sigma TPC vs TOF; n#sigma_{TPC}^{K}; n#sigma_{TOF}^{K}", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});
          rQA.add("tpc_tof_nsigma_pion_mc", "n#sigma TPC vs TOF; n#sigma_{TPC}^{#pi}; n#sigma_{TOF}^{#pi}", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});

          // Rapidity, pseudorapisdity
          rQA.add("hEta_all_after_mc", "Pseudorapidity of all tracks after track selection; #eta; Counts", kTH1F, {{400, -1.0, 1.0}});
          rQA.add("hphi_all_after_mc", "phi of all tracks after track selection; #eta; Counts", kTH1F, {{400, -4.0, 4.0}});

          rQA.add("hEta_ka_mc", "Pseudorapidity of selected Kaons; #eta; Counts", kTH1F, {{400, -1.0, 1.0}});
          rQA.add("hRap_ka_mc", "Rapidity of selected Kaons; y; Counts", kTH1F, {{400, -1.0, 1.0}});

          rQA.add("hEta_pi_mc", "Pseudorapidity of selected Pions; #eta; Counts", kTH1F, {{400, -1.0, 1.0}});
          rQA.add("hRap_pi_mc", "Rapidity of selected Pions; y; Counts", kTH1F, {{400, -1.0, 1.0}});

          // Detector signals
          rQA.add("FT0A_0_mc", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
          rQA.add("FT0A_1_mc", "T0A amplitude", kTH1F, {{20000, 0.0, 20000.0}});
          rQA.add("FT0C_0_mc", "T0C amplitude", kTH1F, {{20000, 0.0, 20000.0}});
          rQA.add("FT0C_1_mc", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
          rQA.add("ZDC_A_0_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
          rQA.add("ZDC_A_1_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
          rQA.add("ZDC_C_0_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
          rQA.add("ZDC_C_1_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
          rQA.add("V0A_0_mc", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
          rQA.add("V0A_1_mc", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
        }
      }
    }
  }

  //_____________________________________________________________________________
  double cosThetaCollinsSoperFrame(ROOT::Math::PxPyPzMVector pair1,
                                   ROOT::Math::PxPyPzMVector pair2,
                                   ROOT::Math::PxPyPzMVector fourpion)
  {
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);

    ROOT::Math::PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    ROOT::Math::PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target

    ROOT::Math::PxPyPzMVector v1 = ROOT::Math::PxPyPzMVector(pair1.Px(), pair1.Py(), pair1.Pz(), pair1.M());
    ROOT::Math::PxPyPzMVector v2 = ROOT::Math::PxPyPzMVector(pair2.Px(), pair2.Py(), pair2.Pz(), pair2.M());
    ROOT::Math::PxPyPzMVector v12 = ROOT::Math::PxPyPzMVector(fourpion.Px(), fourpion.Py(), fourpion.Pz(), fourpion.M());

    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1Cm{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2Cm{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1Cm{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2Cm{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxisCs{((beam1Cm.Unit() - beam2Cm.Unit()).Unit())};

    double cosThetaCs = zaxisCs.Dot((v1Cm));
    return cosThetaCs;
  }

  double phiCollinsSoperFrame(ROOT::Math::PxPyPzMVector pair1, ROOT::Math::PxPyPzMVector pair2, ROOT::Math::PxPyPzMVector fourpion)
  {
    // Half of the energy per pair of the colliding nucleons.
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);

    ROOT::Math::PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    ROOT::Math::PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target

    ROOT::Math::PxPyPzMVector v1 = ROOT::Math::PxPyPzMVector(pair1.Px(), pair1.Py(), pair1.Pz(), pair1.M());
    ROOT::Math::PxPyPzMVector v2 = ROOT::Math::PxPyPzMVector(pair2.Px(), pair2.Py(), pair2.Pz(), pair2.M());
    ROOT::Math::PxPyPzMVector v12 = ROOT::Math::PxPyPzMVector(fourpion.Px(), fourpion.Py(), fourpion.Pz(), fourpion.M());

    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1Cm{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2Cm{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1Cm{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2Cm{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxisCs{((beam1Cm.Unit() - beam2Cm.Unit()).Unit())};
    ROOT::Math::XYZVectorF yaxisCs{(beam1Cm.Cross(beam2Cm)).Unit()};
    ROOT::Math::XYZVectorF xaxisCs{(yaxisCs.Cross(zaxisCs)).Unit()};

    double phi = std::atan2(yaxisCs.Dot(v1Cm), xaxisCs.Dot(v1Cm));
    return phi;
  }

  template <typename C>
  bool isGoodRCTflag(const C& coll)
  {
    switch (cutRCTflag) {
      case 1:
        return sgSelector.isCBTOk(coll);
      case 2:
        return sgSelector.isCBTZdcOk(coll);
      case 3:
        return sgSelector.isCBTHadronOk(coll);
      case 4:
        return sgSelector.isCBTHadronZdcOk(coll);
      default:
        return true;
    }
  }

  template <typename C>
  std::pair<bool, int> selectionEvent(const C& collision, bool fillHist = false)
  {
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 0);

    // Gapside logic
    int gapSide = collision.gapSide();
    float fitCut[5] = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    int truegapSide = sgSelector.trueGap(collision, fitCut[0], fitCut[1], fitCut[2], zdcCut);

    if (fillHist) {
      registry.fill(HIST("GapSide"), gapSide);
      registry.fill(HIST("TrueGapSide"), truegapSide);
    }

    gapSide = truegapSide;

    if (gapSide < SingleGapA || gapSide > DoubleGap)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 1);

    if (std::abs(collision.posZ()) > vzCut)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 2);

    if (useOccCut && (std::abs(collision.occupancyInTime()) > confgOccCut))
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 3);

    if (useHadronicRateCut && (std::abs(collision.hadronicRate()) > confgHadronicRateMax || std::abs(collision.hadronicRate()) < confgHadronicRateMin))
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 4);

    if (useTrs && collision.trs() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 5);

    if (useTrofs && collision.trofs() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 6);

    if (useHmpr && collision.hmpr() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 7);

    if (useTfb && collision.tfb() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 8);

    if (useItsrofb && collision.itsROFb() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 9);

    if (useSbp && collision.sbp() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 10);

    if (useZvtxftovpv && collision.zVtxFT0vPV() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 11);

    if (useVtxItsTpc && collision.vtxITSTPC() != 1)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 12);

    if (!isGoodRCTflag(collision))
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 13);

    if (upcflag != -1 && collision.flags() != upcflag)
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 14);

    if (usenumContrib && (collision.numContrib() < mintrack || collision.numContrib() > maxtrack))
      return {false, gapSide};
    if (fillHist)
      registry.fill(HIST("hEventCutFlow"), 15);

    return {true, gapSide};
  }

  template <typename T>
  bool selectionPIDKaon1(const T& candidate)
  {
    auto pt = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py());
    //  float pt1, pt2, pt3 , nsigmatpc_cut1, nsigmatpc_cut2, nsigmatpc_cut3;
    if (useTof && pt < pt1 && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut1) {
      return true;
    }
    if (useTof && pt >= pt1 && pt < pt2 && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut2) {
      return true;
    }
    if (useTof && pt >= pt2 && pt < pt3 && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut3) {
      return true;
    }
    if (useTof && pt >= pt3 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && ccut && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < nsigmaTofCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && !ccut && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut && std::abs(candidate.tofNSigmaKa()) < nsigmaTofCut1) {
      return true;
    }

    if (!useTof && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPIDPion1(const T& candidate)
  {
    auto pt = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py());

    if (useTof && pt < pt1 && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut1) {
      return true;
    }
    if (useTof && pt >= pt1 && pt < pt2 && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut2) {
      return true;
    }
    if (useTof && pt >= pt2 && pt < pt3 && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut3) {
      return true;
    }
    if (useTof && pt >= pt3 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && ccut && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < nsigmaTofCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && !ccut && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut && std::abs(candidate.tofNSigmaPi()) < nsigmaTofCut1) {
      return true;
    }
    if (!useTof && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut) {
      return true;
    }
    return false;
  }

  //------------------------------------------------------------------------------------------------------

  using UDtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& collision, UDtracksfull const& tracks)
  {
    if (qa)
      rQA.fill(HIST("hOcc_before"), collision.occupancyInTime());

    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;
    ROOT::Math::PxPyPzMVector v01;

    ROOT::Math::PxPyPzMVector phiv;
    ROOT::Math::PxPyPzMVector phiv1;

    std::vector<ROOT::Math::PxPyPzMVector> onlyPionTracksp;
    std::vector<decltype(tracks.begin())> rawPionTracksp;

    std::vector<ROOT::Math::PxPyPzMVector> onlyPionTrackspm;
    std::vector<decltype(tracks.begin())> rawPionTrackspm;

    std::vector<ROOT::Math::PxPyPzMVector> onlyPionTracksn;
    std::vector<decltype(tracks.begin())> rawPionTracksn;

    std::vector<float> parameters = {pvCut, dcazCut, dcaxyCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};

    auto [eventSelected, gapSide] = selectionEvent(collision, true);
    if (!eventSelected)
      return;

    int mult = collision.numContrib();

    if (gapSide == SingleGapA) {
      registry.fill(HIST("gap_mult0"), mult);
    }
    if (gapSide == SingleGapC) {
      registry.fill(HIST("gap_mult1"), mult);
    }
    if (gapSide == DoubleGap) {
      registry.fill(HIST("gap_mult2"), mult);
    }

    if (qa)
      rQA.fill(HIST("hOcc_after"), collision.occupancyInTime());

    int mult0 = 0;
    int mult1 = 0;
    int mult2 = 0;
    if (qa) {
      rQA.fill(HIST("hVertexX"), collision.posX());
      rQA.fill(HIST("hVertexY"), collision.posY());
      rQA.fill(HIST("hVertexZ"), collision.posZ());
    }

    for (const auto& track1 : tracks) {
      if (qa) {
        rQA.fill(HIST("hDcaxy_all_before"), track1.dcaXY());
        rQA.fill(HIST("hDcaz_all_before"), track1.dcaZ());
      }

      if (!trackselector(track1, parameters))
        continue;

      v0.SetCoordinates(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassPionCharged);

      if (qa) {
        rQA.fill(HIST("hDcaxy_all_after"), track1.dcaXY());
        rQA.fill(HIST("hDcaz_all_after"), track1.dcaZ());
        rQA.fill(HIST("hEta_all_after"), v0.Eta());
      }

      if (selectionPIDPion1(track1)) {
        onlyPionTrackspm.push_back(v0);
        rawPionTrackspm.push_back(track1);
        if (track1.sign() == 1) {
          onlyPionTracksp.push_back(v0);
          rawPionTracksp.push_back(track1);
        }
        if (track1.sign() == -1) {
          onlyPionTracksn.push_back(v0);
          rawPionTracksn.push_back(track1);
        }
      }
      if (gapSide == SingleGapA) {
        mult0++;
      }
      if (gapSide == SingleGapC) {
        mult1++;
      }
      if (gapSide == DoubleGap) {
        mult2++;
      }

      if (qa) {
        rQA.fill(HIST("tpc_dedx"), v0.P(), track1.tpcSignal());
        rQA.fill(HIST("tof_beta"), v0.P(), track1.beta());
        rQA.fill(HIST("tof_nsigma_kaon_all"), v0.Pt(), track1.tofNSigmaKa());
        rQA.fill(HIST("tof_nsigma_pion_all"), v0.Pt(), track1.tofNSigmaPi());
        rQA.fill(HIST("tpc_nsigma_kaon_all"), v0.Pt(), track1.tpcNSigmaKa());
        rQA.fill(HIST("tpc_nsigma_pion_all"), v0.Pt(), track1.tpcNSigmaPi());

        if (selectionPIDKaon1(track1)) {
          rQA.fill(HIST("tpc_dedx_kaon"), v0.P(), track1.tpcSignal());
          rQA.fill(HIST("tof_beta_kaon"), v0.P(), track1.beta());
          rQA.fill(HIST("tpc_nsigma_kaon"), v0.Pt(), track1.tpcNSigmaKa());
          rQA.fill(HIST("tof_nsigma_kaon"), v0.Pt(), track1.tofNSigmaKa());
          rQA.fill(HIST("tpc_tof_nsigma_kaon"), track1.tpcNSigmaKa(), track1.tofNSigmaKa());
          rQA.fill(HIST("hEta_ka"), v0.Eta());
          rQA.fill(HIST("hRap_ka"), v0.Rapidity());
          rQA.fill(HIST("hDcaxy_ka"), track1.dcaXY());
          rQA.fill(HIST("hDcaz_ka"), track1.dcaZ());
        }

        if (selectionPIDPion1(track1)) {
          rQA.fill(HIST("tpc_dedx_pion"), v0.P(), track1.tpcSignal());
          rQA.fill(HIST("tof_beta_pion"), v0.P(), track1.beta());
          rQA.fill(HIST("tpc_nsigma_pion"), v0.Pt(), track1.tpcNSigmaPi());
          rQA.fill(HIST("tof_nsigma_pion"), v0.Pt(), track1.tofNSigmaPi());
          rQA.fill(HIST("tpc_tof_nsigma_pion"), track1.tpcNSigmaPi(), track1.tofNSigmaPi());
          rQA.fill(HIST("hEta_pi"), v0.Eta());
          rQA.fill(HIST("hRap_pi"), v0.Rapidity());
          rQA.fill(HIST("hDcaxy_pi"), track1.dcaXY());
          rQA.fill(HIST("hDcaz_pi"), track1.dcaZ());
        }
      }
    }
    if (gapSide == SingleGapA) {
      if (useMultCut && (mult0 < mintrack || mult0 > maxtrack))
        return;
      registry.fill(HIST("mult_0"), mult0);
    }
    if (gapSide == SingleGapC) {
      if (useMultCut && (mult1 < mintrack || mult1 > maxtrack))
        return;
      registry.fill(HIST("mult_1"), mult1);
    }
    if (qa) {
      if (gapSide == SingleGapA) {
        rQA.fill(HIST("V0A_0"), collision.totalFV0AmplitudeA());
        rQA.fill(HIST("FT0A_0"), collision.totalFT0AmplitudeA());
        rQA.fill(HIST("FT0C_0"), collision.totalFT0AmplitudeC());
        rQA.fill(HIST("ZDC_A_0"), collision.energyCommonZNA());
        rQA.fill(HIST("ZDC_C_0"), collision.energyCommonZNC());
      }
      if (gapSide == SingleGapC) {
        rQA.fill(HIST("V0A_1"), collision.totalFV0AmplitudeA());
        rQA.fill(HIST("FT0A_1"), collision.totalFT0AmplitudeA());
        rQA.fill(HIST("FT0C_1"), collision.totalFT0AmplitudeC());
        rQA.fill(HIST("ZDC_A_1"), collision.energyCommonZNA());
        rQA.fill(HIST("ZDC_C_1"), collision.energyCommonZNC());
      }
    }
    for (const auto& [t0, t1] : combinations(tracks, tracks)) {
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;

      if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
        // Apply kaon hypothesis and create pairs
        v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
        v01 = v0 + v1;
        if (rapiditycut && std::abs(v01.Rapidity()) > rapiditycutvalue)
          continue;

        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == SingleGapA) {
            registry.fill(HIST("os_KK_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == SingleGapC) {
            registry.fill(HIST("os_KK_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
            registry.fill(HIST("os_KK_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }
        // samesignpair
        if (t0.sign() == t1.sign()) {
          if (gapSide == SingleGapA) {
            if (t0.sign() < 0)
              registry.fill(HIST("os_KK_lsMM_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
            else
              registry.fill(HIST("os_KK_lsPP_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == SingleGapC) {
            if (t0.sign() < 0)
              registry.fill(HIST("os_KK_lsMM_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
            else
              registry.fill(HIST("os_KK_lsPP_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
            registry.fill(HIST("os_KK_ls_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }

        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            registry.fill(HIST("hRotation"), rotangle);

            auto rotkaonPx = t0.px() * std::cos(rotangle) - t0.py() * std::sin(rotangle);
            auto rotkaonPy = t0.px() * std::sin(rotangle) + t0.py() * std::cos(rotangle);

            v0.SetCoordinates(rotkaonPx, rotkaonPy, t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
            v01 = v0 + v1;
            if (rapiditycut && std::abs(v01.Rapidity()) > rapiditycutvalue)
              continue;

            if (t0.sign() != t1.sign()) {
              if (gapSide == SingleGapA) {
                registry.fill(HIST("os_KK_rot_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == SingleGapC) {
                registry.fill(HIST("os_KK_rot_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
                registry.fill(HIST("os_KK_rot_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
    }
    for (const auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;
      if (t0.globalIndex() == t1.globalIndex())
        continue;
      if (rho && selectionPIDProton(t0, useTof, nsigmaTpcCut, nsigmaTofCut) && selectionPIDKaon1(t1)) {
        v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassProton);
        v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
        v01 = v0 + v1;
        if (rapiditycut && std::abs(v01.Rapidity()) > rapiditycutvalue)
          continue;

        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == SingleGapA) {
            registry.fill(HIST("os_pp_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == SingleGapC) {
            registry.fill(HIST("os_pp_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
            registry.fill(HIST("os_pp_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        } // same sign pair
        if (t0.sign() == t1.sign()) {
          if (gapSide == SingleGapA) {
            registry.fill(HIST("os_pp_ls_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == SingleGapC) {
            registry.fill(HIST("os_pp_ls_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
            registry.fill(HIST("os_pp_ls_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }
      }
      if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
        if (kaoncut && t0.tpcNSigmaPi() < pionNsigmaCut)
          continue;
        v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
        v01 = v0 + v1;
        if (rapiditycut && std::abs(v01.Rapidity()) > rapiditycutvalue)
          continue;

        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == SingleGapA) {
            registry.fill(HIST("os_pk_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == SingleGapC) {
            registry.fill(HIST("os_pk_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
            registry.fill(HIST("os_pk_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        } // same sign pair
        if (t0.sign() == t1.sign()) {
          if (gapSide == SingleGapA) {
            if (t0.sign() < 0)
              registry.fill(HIST("os_pk_lsMM_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
            else
              registry.fill(HIST("os_pk_lsPP_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == SingleGapC) {
            if (t0.sign() < 0)
              registry.fill(HIST("os_pk_lsMM_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
            else
              registry.fill(HIST("os_pk_lsPP_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
            registry.fill(HIST("os_pk_ls_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            registry.fill(HIST("hRotation"), rotangle);

            auto rotkaonPx = t0.px() * std::cos(rotangle) - t0.py() * std::sin(rotangle);
            auto rotkaonPy = t0.px() * std::sin(rotangle) + t0.py() * std::cos(rotangle);

            v0.SetCoordinates(rotkaonPx, rotkaonPy, t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
            v01 = v0 + v1;
            if (rapiditycut && std::abs(v01.Rapidity()) > rapiditycutvalue)
              continue;

            if (t0.sign() != t1.sign()) {
              if (gapSide == SingleGapA) {
                registry.fill(HIST("os_pk_rot_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == SingleGapC) {
                registry.fill(HIST("os_pk_rot_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (exclusive && gapSide == DoubleGap && mult2 == numTwoTracks) {
                registry.fill(HIST("os_pk_rot_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
    }
    if (fourpion) {
      if (gapSide == DoubleGap && mult2 == numFourTracks) {
        ROOT::Math::PxPyPzMVector pair1, pair2, pair3, pair4;
        if (static_cast<int>(onlyPionTracksp.size()) == numTwoTracks && static_cast<int>(onlyPionTracksn.size()) == numTwoTracks) {
          ROOT::Math::PxPyPzMVector k1 = onlyPionTracksp.at(0);
          ROOT::Math::PxPyPzMVector k2 = onlyPionTracksp.at(1);
          ROOT::Math::PxPyPzMVector k3 = onlyPionTracksn.at(0);
          ROOT::Math::PxPyPzMVector k4 = onlyPionTracksn.at(1);
          phiv = k1 + k2 + k3 + k4;
          pair1 = k1 + k3;
          pair2 = k2 + k4;
          pair3 = k1 + k4;
          pair4 = k2 + k3;
          registry.fill(HIST("os_pppp_pT_2"), phiv.M(), phiv.Pt(), phiv.Rapidity());
          registry.fill(HIST("os_pp_vs_pp_mass"), pair1.M(), pair2.M());
          registry.fill(HIST("os_pp_vs_pp_pt"), pair1.Pt(), pair2.Pt());
          auto costhetaPair = cosThetaCollinsSoperFrame(pair1, pair2, phiv);
          auto phiPair = 1. * o2::constants::math::PI + phiCollinsSoperFrame(pair1, pair2, phiv);
          registry.fill(HIST("phi_dis"), phiPair);
          registry.fill(HIST("costheta_dis"), costhetaPair);
          registry.fill(HIST("costheta_vs_phi"), costhetaPair, phiPair);
          registry.fill(HIST("os_pp_vs_pp_mass1"), pair3.M(), pair4.M());
          registry.fill(HIST("os_pp_vs_pp_pt1"), pair3.Pt(), pair4.Pt());
          auto costhetaPair1 = cosThetaCollinsSoperFrame(pair3, pair4, phiv);
          auto phiPair1 = 1. * o2::constants::math::PI + phiCollinsSoperFrame(pair3, pair4, phiv);
          registry.fill(HIST("phi_dis1"), phiPair1);
          registry.fill(HIST("costheta_dis1"), costhetaPair1);
          registry.fill(HIST("costheta_vs_phi1"), costhetaPair1, phiPair1);
        }
        if (static_cast<int>(onlyPionTracksp.size()) != numTwoTracks && static_cast<int>(onlyPionTracksn.size()) != numTwoTracks) {
          if (static_cast<int>(onlyPionTracksp.size()) + static_cast<int>(onlyPionTracksn.size()) != numFourTracks)
            return;
          ROOT::Math::PxPyPzMVector l1 = onlyPionTrackspm.at(0);
          ROOT::Math::PxPyPzMVector l2 = onlyPionTrackspm.at(1);
          ROOT::Math::PxPyPzMVector l3 = onlyPionTrackspm.at(2);
          ROOT::Math::PxPyPzMVector l4 = onlyPionTrackspm.at(3);
          phiv1 = l1 + l2 + l3 + l4;
          registry.fill(HIST("os_pppp_pT_2_ls"), phiv1.M(), phiv1.Pt(), phiv1.Rapidity());
        }
      }
    }
  }
  PROCESS_SWITCH(SginclusivePhiKstarSD, process, "Process unlike event", false);

  using UDCollisionsFull1 = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  SliceCache cache;
  Partition<UDtracksfull> posTracks = aod::udtrack::sign > 0;
  Partition<UDtracksfull> negTracks = aod::udtrack::sign < 0;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for bin"};
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  void mixprocess(UDCollisionsFull1 const& collisions, UDtracksfull const& /*track*/)
  {
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;
    ROOT::Math::PxPyPzMVector v01;

    std::vector<float> parameters = {pvCut, dcazCut, dcaxyCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};

    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass}, true};

    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {

      auto [eventSelected1, gapSide1] = selectionEvent(collision1, false);
      auto [eventSelected2, gapSide2] = selectionEvent(collision2, false);
      if (!eventSelected1 || !eventSelected2)
        continue;

      if (gapSide1 != gapSide2)
        continue;

      auto posThisColl = posTracks->sliceByCached(aod::udtrack::udCollisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::udtrack::udCollisionId, collision2.globalIndex(), cache);
      //      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
      for (const auto& [track1, track2] : o2::soa::combinations(posThisColl, negThisColl)) {
        if (!trackselector(track1, parameters) || !trackselector(track2, parameters))
          continue;
        if (phi && selectionPIDKaon1(track1) && selectionPIDKaon1(track2)) {
          v0.SetCoordinates(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassKaonCharged);
          v1.SetCoordinates(track2.px(), track2.py(), track2.pz(), o2::constants::physics::MassKaonCharged);
          v01 = v0 + v1;
          if (rapiditycut && std::abs(v01.Rapidity()) > rapiditycutvalue)
            continue;
          // Opposite sign pairs
          if (track1.sign() != track2.sign()) {
            if (gapSide1 == SingleGapA) {
              registry.fill(HIST("os_KK_mix_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
            }
            if (gapSide1 == SingleGapC) {
              registry.fill(HIST("os_KK_mix_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
            }
          }
        }
      }
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
        if (!trackselector(track1, parameters) || !trackselector(track2, parameters))
          continue;
        if (track1.globalIndex() == track2.globalIndex())
          continue;
        if (kstar && selectionPIDKaon1(track1) && selectionPIDPion1(track2)) {
          v0.SetCoordinates(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassKaonCharged);
          v1.SetCoordinates(track2.px(), track2.py(), track2.pz(), o2::constants::physics::MassPionCharged);
          v01 = v0 + v1;
          if (rapiditycut && std::abs(v01.Rapidity()) > rapiditycutvalue)
            continue;
          // Opposite sign pairs
          if (track1.sign() != track2.sign()) {
            if (gapSide1 == SingleGapA) {
              registry.fill(HIST("os_pk_mix_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
            }
            if (gapSide1 == SingleGapC) {
              registry.fill(HIST("os_pk_mix_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(SginclusivePhiKstarSD, mixprocess, "Process Mixed event", false);

  // define abbreviations  , aod::UDCollisions_001,
  using CCs = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDMcCollsLabels, aod::UDZdcsReduced>;
  using CC = CCs::iterator;
  // using TCs = soa::Join<aod::UDTracks>;
  using TCs = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksFlags, aod::UDTracksPID, aod::UDMcTrackLabels>;
  using TC = TCs::iterator;

  PresliceUnsorted<aod::UDMcParticles> partPerMcCollision = aod::udmcparticle::udMcCollisionId;
  PresliceUnsorted<CCs> colPerMcCollision = aod::udcollision::udMcCollisionId;
  PresliceUnsorted<TCs> trackPerMcParticle = aod::udmctracklabel::udMcParticleId;

  void processMCTruth(aod::UDMcCollisions const& mccollisions, CCs const& collisions, aod::UDMcParticles const& McParts, TCs const& /*tracks*/)
  {
    // number of McCollisions in DF
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;
    ROOT::Math::PxPyPzMVector v01;
    ROOT::Math::PxPyPzMVector vkstar;
    ROOT::Math::PxPyPzMVector vphi;
    for (const auto& mccollision : mccollisions) {
      if (mccollision.generatorsID() != generatedId)
        continue;
      registry.get<TH1>(HIST("MC/Stat"))->Fill(0., 1.);
      // get reconstructed collision which belongs to mccollision
      auto colSlice = collisions.sliceBy(colPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/recCols"))->Fill(colSlice.size(), 1.);
      if (reconstruction && colSlice.size() < 1)
        continue;
      // get McParticles which belong to mccollision
      auto partSlice = McParts.sliceBy(partPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/nParts"))->Fill(partSlice.size(), 1.);
      for (const auto& trk1 : partSlice) {
        if (std::abs(trk1.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
          vkstar.SetCoordinates(trk1.px(), trk1.py(), trk1.pz(), o2::constants::physics::MassK0Star892);
          auto kDaughters = trk1.daughters_as<aod::UDMcParticles>();
          const size_t kExpectedNumberOfDaughters = 2;
          if (kDaughters.size() != kExpectedNumberOfDaughters) {
            continue;
          }
          auto daughtp = false;
          auto daughtm = false;
          for (const auto& kCurrentDaughter : kDaughters) {
            if (!kCurrentDaughter.isPhysicalPrimary()) {
              continue;
            }
            if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kKPlus) {
              daughtp = true;
              v0.SetCoordinates(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassKaonCharged);
            } else if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kPiPlus) {
              daughtm = true;
              v1.SetCoordinates(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassPionCharged);
            }
          }
          if (daughtp && daughtm) {
            v01 = v0 + v1;
            if (std::abs(v01.Rapidity()) < 0.5) {
              registry.get<TH3>(HIST("MC/accMPtRap_kstar_G"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
            }
            if (std::abs(vkstar.Rapidity()) < 0.5) {
              //  registry.get<TH3>(HIST("MC/accMPtRap_kstar_G1"))->Fill(vkstar.M(), vkstar.Pt(), vkstar.Rapidity(), 1.);
            }
          }
        }
        // phi
        if (std::abs(trk1.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
          if (trk1.has_mothers()) {
            continue;
          }
          vphi.SetCoordinates(trk1.px(), trk1.py(), trk1.pz(), o2::constants::physics::MassPhi);
          auto kDaughters = trk1.daughters_as<aod::UDMcParticles>();
          const size_t kExpectedNumberOfDaughters1 = 2;
          if (kDaughters.size() != kExpectedNumberOfDaughters1) {
            continue;
          }
          auto daughtp1 = false;
          auto daughtm1 = false;
          for (const auto& kCurrentDaughter : kDaughters) {
            if (!kCurrentDaughter.isPhysicalPrimary()) {
              continue;
            }
            if (kCurrentDaughter.pdgCode() == PDG_t::kKPlus) {
              daughtp1 = true;
              v0.SetCoordinates(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassKaonCharged);
            } else if (kCurrentDaughter.pdgCode() == PDG_t::kKMinus) {
              daughtm1 = true;
              v1.SetCoordinates(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassKaonCharged);
            }
          }
          if (daughtp1 && daughtm1) {
            v01 = v0 + v1;
            if (std::abs(v01.Rapidity()) < 0.5) {
              registry.get<TH3>(HIST("MC/accMPtRap_phi_G"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
            }
            if (std::abs(vphi.Rapidity()) < 0.5) {
              registry.get<TH3>(HIST("MC/accMPtRap_phi_G2"))->Fill(vphi.M(), vphi.Pt(), vphi.Rapidity(), 1.);
            }
          }
        }
      }
      /*
      for (const auto& [tr1, tr2] : combinations(partSlice, partSlice)) {
        if ((tr1.pdgCode() == kKPlus && tr2.pdgCode() == kPiMinus) || (tr1.pdgCode() == kKMinus && tr2.pdgCode() == kPiPlus) || (tr1.pdgCode() == kPiPlus && tr2.pdgCode() == kKMinus) || (tr1.pdgCode() == kPiMinus && tr2.pdgCode() == kKPlus)) {
          if (std::abs(tr1.pdgCode()) == kKPlus) {
            v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassPionCharged);
          } else {
            v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassPionCharged);
            v1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
          }
          if (!tr1.isPhysicalPrimary() || !tr2.isPhysicalPrimary())
            continue;
          v01 = v0 + v1;
          if (tr1.globalIndex() + 1 != tr2.globalIndex()) {
            registry.get<TH1>(HIST("MC/genM_1_k"))->Fill(v01.M(), 1.);
          }
          if (std::abs(tr1.globalIndex() - tr2.globalIndex()) != 1)
            continue;
          bool flag = false;
          bool flag1 = false;
          if (tr1.has_mothers() && tr2.has_mothers()) {
            for (const auto& mother : tr1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
                vkstar.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassK0Star892);
    // registry.get<TH3>(HIST("MC/accMPtRap_kstar_G"))->Fill(vkstar.M(), vkstar.Pt(), vkstar.Rapidity(), 1.);
                flag = true;
              }
            }
            for (const auto& mother1 : tr2.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
                flag1 = true;
              }
            }
          }
          if (flag && flag1) {
            registry.get<TH2>(HIST("MC/genMPt_k"))->Fill(v01.M(), v01.Pt(), 1.);
          }
          registry.get<TH1>(HIST("MC/genRap_k"))->Fill(v01.Rapidity(), 1.);
          registry.get<TH1>(HIST("MC/genM_k"))->Fill(v01.M(), 1.);
          registry.get<TH1>(HIST("MC/accRap_k"))->Fill(v01.Rapidity(), 1.);
          registry.get<TH2>(HIST("MC/accMPt_k"))->Fill(v01.M(), v01.Pt(), 1.);
          registry.get<TH1>(HIST("MC/accM_k"))->Fill(v01.M(), 1.);
          registry.get<TH3>(HIST("MC/accMPtRap_k"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
        }

        if (std::abs(tr1.pdgCode()) != kKPlus || std::abs(tr2.pdgCode()) != kKPlus)
          continue;
        v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
        if (tr1.pdgCode() == tr2.pdgCode())
          continue;
        v01 = v0 + v1;
        if (!tr1.isPhysicalPrimary() || !tr2.isPhysicalPrimary())
          continue;
        if (tr1.globalIndex() + 1 != tr2.globalIndex()) {
          registry.get<TH1>(HIST("MC/genM_1"))->Fill(v01.M(), 1.);
        }
        if (std::abs(tr1.globalIndex() - tr2.globalIndex()) != 1)
          continue;
        bool flag = false;
        bool flag1 = false;
        if (tr1.has_mothers() && tr2.has_mothers()) {
          for (const auto& mother : tr1.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
              vphi.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassPhi);
        //  registry.get<TH3>(HIST("MC/accMPtRap_phi_G"))->Fill(vphi.M(), vphi.Pt(), vphi.Rapidity(), 1.);
              flag = true;
            }
          }
          for (const auto& mother1 : tr2.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
              flag1 = true;
            }
          }
        }
        if (flag && flag1) {
          registry.get<TH2>(HIST("MC/genMPt"))->Fill(v01.M(), v01.Pt(), 1.);
        }
        // registry.get<TH1>(HIST("MC/genRap"))->Fill(v01.Rapidity(), 1.);
        // registry.get<TH1>(HIST("MC/genM"))->Fill(v01.M(), 1.);
        registry.get<TH1>(HIST("MC/accRap"))->Fill(v01.Rapidity(), 1.);
        registry.get<TH2>(HIST("MC/accMPt"))->Fill(v01.M(), v01.Pt(), 1.);
        registry.get<TH1>(HIST("MC/accM"))->Fill(v01.M(), 1.);
        registry.get<TH3>(HIST("MC/accMPtRap"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
      }
      // compute the difference between generated and reconstructed particle momentum
      for (const auto& McPart : partSlice) {
        auto trackSlice = tracks.sliceBy(trackPerMcParticle, McPart.globalIndex());
        registry.get<TH1>(HIST("MC/nRecTracks"))->Fill(trackSlice.size(), 1.);
        // compute momentum difference between MCTruth and Reconstruction
        if (trackSlice.size() > 0) {
          for (const auto& track : trackSlice) {
            // std::cout<<trackSlice.size()<<std::endl;
            auto pTrack = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
            auto pPart = std::sqrt(McPart.px() * McPart.px() + McPart.py() * McPart.py() + McPart.pz() * McPart.pz());
            auto pDiff = pTrack - pPart;
            registry.get<TH2>(HIST("MC/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
          }
        } else {
          registry.get<TH2>(HIST("MC/pDiff"))->Fill(-5.9, -1, 1.);
        }
      }
      */
    }
  }
  PROCESS_SWITCH(SginclusivePhiKstarSD, processMCTruth, "Process MC truth", true);
  // ...............................................................................................................
  void processReco(CC const& collision, TCs const& tracks, aod::UDMcCollisions const& /*mccollisions*/, aod::UDMcParticles const& McParts)
  {
    registry.fill(HIST("Reco/hEventCutFlowMC"), 0);
    ROOT::Math::PxPyPzMVector vkstar;
    ROOT::Math::PxPyPzMVector vphi;

    if (!collision.has_udMcCollision())
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 1);

    if (qaMC)
      rQA.fill(HIST("hOcc_before_mc"), collision.occupancyInTime());

    auto mccoll = collision.udMcCollision();
    if (mccoll.generatorsID() != generatedId)
      return;
    if (upcflag != -1 && collision.flags() == upcflag) {
      auto partSlice = McParts.sliceBy(partPerMcCollision, mccoll.globalIndex());
      for (const auto& trk1 : partSlice) {
        if (std::abs(trk1.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
          vkstar.SetCoordinates(trk1.px(), trk1.py(), trk1.pz(), o2::constants::physics::MassK0Star892);
          if (std::abs(vkstar.Rapidity()) < 0.5) {
            registry.get<TH3>(HIST("MC/accMPtRap_kstar_G1"))->Fill(vkstar.M(), vkstar.Pt(), vkstar.Rapidity(), 1.);
          }
        }
        // phi
        if (std::abs(trk1.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
          if (trk1.has_mothers()) {
            continue;
          }
          vphi.SetCoordinates(trk1.px(), trk1.py(), trk1.pz(), o2::constants::physics::MassPhi);
          if (std::abs(vphi.Rapidity()) < 0.5) {
            registry.get<TH3>(HIST("MC/accMPtRap_phi_G1"))->Fill(vphi.M(), vphi.Pt(), vphi.Rapidity(), 1.);
          }
        }
      }
    }

    registry.fill(HIST("Reco/hEventCutFlowMC"), 2);

    float fitCut[5] = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    std::vector<float> parameters = {pvCut, dcazCut, dcaxyCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};
    int truegapSide = sgSelector.trueGap(collision, fitCut[0], fitCut[1], fitCut[2], zdcCut);
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(4.0, 1.);
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(truegapSide, 1.);
    if (truegapSide != gapsideMC)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 3);
    //    Partition<TCs> pvContributors = aod::udtrack::isPVContributor == true;
    // pvContributors.bindTable(tracks);

    if (std::abs(collision.posZ()) > vzCut)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 4);

    if (useOccCut && (std::abs(collision.occupancyInTime()) > confgOccCut))
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 5);

    if (useHadronicRateCut && (std::abs(collision.hadronicRate()) > confgHadronicRateMax || std::abs(collision.hadronicRate()) < confgHadronicRateMin))
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 6);

    if (useTrs && collision.trs() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 7);

    if (useTrofs && collision.trofs() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 8);

    if (useHmpr && collision.hmpr() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 9);

    if (useTfb && collision.tfb() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 10);

    if (useItsrofb && collision.itsROFb() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 11);

    if (useSbp && collision.sbp() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 12);

    if (useZvtxftovpv && collision.zVtxFT0vPV() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 13);

    if (useVtxItsTpc && collision.vtxITSTPC() != 1)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 14);

    if (!isGoodRCTflag(collision))
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 15);

    if (upcflag != -1 && collision.flags() != upcflag)
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 16);

    if (usenumContrib && (collision.numContrib() < mintrack || collision.numContrib() > maxtrack))
      return;
    registry.fill(HIST("Reco/hEventCutFlowMC"), 17);

    if (qaMC)
      rQA.fill(HIST("hOcc_after_mc"), collision.occupancyInTime());

    // registry.get<TH1>(HIST("Reco/nPVContributors"))->Fill(pvContributors.size(), 1.);
    // ROOT::Math::PxPyPzMVector vphi;
    // ROOT::Math::PxPyPzMVector vkstar;
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector vr0;
    ROOT::Math::PxPyPzMVector vr1;
    ROOT::Math::PxPyPzMVector vr01;
    ROOT::Math::PxPyPzMVector vr0g;
    ROOT::Math::PxPyPzMVector vr1g;
    ROOT::Math::PxPyPzMVector vr01g;
    int t1 = 0;
    if (truegapSide == SingleGapA) {
      if (qaMC) {
        rQA.fill(HIST("V0A_0_mc"), collision.totalFV0AmplitudeA());
        rQA.fill(HIST("FT0A_0_mc"), collision.totalFT0AmplitudeA());
        rQA.fill(HIST("FT0C_0_mc"), collision.totalFT0AmplitudeC());
        rQA.fill(HIST("ZDC_A_0_mc"), collision.energyCommonZNA());
        rQA.fill(HIST("ZDC_C_0_mc"), collision.energyCommonZNC());
      }
    }
    if (truegapSide == SingleGapC) {
      if (qaMC) {
        rQA.fill(HIST("V0A_1_mc"), collision.totalFV0AmplitudeA());
        rQA.fill(HIST("FT0A_1_mc"), collision.totalFT0AmplitudeA());
        rQA.fill(HIST("FT0C_1_mc"), collision.totalFT0AmplitudeC());
        rQA.fill(HIST("ZDC_A_1_mc"), collision.energyCommonZNA());
        rQA.fill(HIST("ZDC_C_1_mc"), collision.energyCommonZNC());
      }
    }
    for (const auto& tr1 : tracks) {
      if (!tr1.has_udMcParticle())
        continue;
      auto mcPart1 = tr1.udMcParticle();

      if (qaMC) {
        rQA.fill(HIST("hDcaxy_all_before_mc"), tr1.dcaXY());
        rQA.fill(HIST("hDcaz_all_before_mc"), tr1.dcaZ());
      }

      registry.get<TH1>(HIST("Reco/tr_dcaz_1"))->Fill(tr1.dcaZ(), 1.);
      registry.get<TH1>(HIST("Reco/tr_dcaxy_1"))->Fill(tr1.dcaXY(), 1.);
      registry.get<TH1>(HIST("Reco/tr_chi2ncl_1"))->Fill(tr1.tpcChi2NCl(), 1.);
      registry.get<TH1>(HIST("Reco/tr_tpcnclfind_1"))->Fill(tr1.tpcNClsFindable(), 1.);
      registry.get<TH1>(HIST("Reco/tr_itsChi2NCl_1"))->Fill(tr1.itsChi2NCl(), 1.);

      if (!trackselector(tr1, parameters))
        continue;

      registry.get<TH1>(HIST("Reco/tr_dcaz_2"))->Fill(tr1.dcaZ(), 1.);
      registry.get<TH1>(HIST("Reco/tr_dcaxy_2"))->Fill(tr1.dcaXY(), 1.);
      registry.get<TH1>(HIST("Reco/tr_chi2ncl_2"))->Fill(tr1.tpcChi2NCl(), 1.);
      registry.get<TH1>(HIST("Reco/tr_tpcnclfind_2"))->Fill(tr1.tpcNClsFindable(), 1.);
      registry.get<TH1>(HIST("Reco/tr_itsChi2NCl_2"))->Fill(tr1.itsChi2NCl(), 1.);
      v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassPionCharged);

      if (qaMC) {
        rQA.fill(HIST("hDcaxy_all_after_mc"), tr1.dcaXY());
        rQA.fill(HIST("hDcaz_all_after_mc"), tr1.dcaZ());
        rQA.fill(HIST("hEta_all_after_mc"), v0.Eta());
        rQA.fill(HIST("hphi_all_after_mc"), v0.Phi());

        rQA.fill(HIST("tpc_dedx_mc"), v0.P(), tr1.tpcSignal());
        rQA.fill(HIST("tof_beta_mc"), v0.P(), tr1.beta());
        rQA.fill(HIST("tof_nsigma_kaon_all_mc"), v0.Pt(), tr1.tofNSigmaKa());
        rQA.fill(HIST("tof_nsigma_pion_all_mc"), v0.Pt(), tr1.tofNSigmaPi());
        rQA.fill(HIST("tpc_nsigma_kaon_all_mc"), v0.Pt(), tr1.tpcNSigmaKa());
        rQA.fill(HIST("tpc_nsigma_pion_all_mc"), v0.Pt(), tr1.tpcNSigmaPi());
      }
      if (selectionPIDKaon1(tr1)) {
        if (qaMC) {
          rQA.fill(HIST("tpc_dedx_kaon_mc"), v0.P(), tr1.tpcSignal());
          rQA.fill(HIST("tof_beta_kaon_mc"), v0.P(), tr1.beta());
          rQA.fill(HIST("tpc_nsigma_kaon_mc"), v0.Pt(), tr1.tpcNSigmaKa());
          rQA.fill(HIST("tof_nsigma_kaon_mc"), v0.Pt(), tr1.tofNSigmaKa());
          rQA.fill(HIST("tpc_tof_nsigma_kaon_mc"), tr1.tpcNSigmaKa(), tr1.tofNSigmaKa());
          rQA.fill(HIST("hEta_ka_mc"), v0.Eta());
          rQA.fill(HIST("hRap_ka_mc"), v0.Rapidity());
          rQA.fill(HIST("hDcaxy_ka_mc"), tr1.dcaXY());
          rQA.fill(HIST("hDcaz_ka_mc"), tr1.dcaZ());
        }
      }
      if (selectionPIDPion1(tr1)) {
        if (qaMC) {
          rQA.fill(HIST("tpc_dedx_pion_mc"), v0.P(), tr1.tpcSignal());
          rQA.fill(HIST("tof_beta_pion_mc"), v0.P(), tr1.beta());
          rQA.fill(HIST("tpc_nsigma_pion_mc"), v0.Pt(), tr1.tpcNSigmaPi());
          rQA.fill(HIST("tof_nsigma_pion_mc"), v0.Pt(), tr1.tofNSigmaPi());
          rQA.fill(HIST("tpc_tof_nsigma_pion_mc"), tr1.tpcNSigmaPi(), tr1.tofNSigmaPi());
          rQA.fill(HIST("hEta_pi_mc"), v0.Eta());
          rQA.fill(HIST("hRap_pi_mc"), v0.Rapidity());
          rQA.fill(HIST("hDcaxy_pi_mc"), tr1.dcaXY());
          rQA.fill(HIST("hDcaz_pi_mc"), tr1.dcaZ());
        }
      }

      t1++;
      vr0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassKaonCharged);
      registry.get<TH1>(HIST("Reco/trpt"))->Fill(vr0.Pt(), 1.);
      registry.get<TH1>(HIST("Reco/treta_k"))->Fill(vr0.Eta(), 1.);
      if (!selectionPIDKaon1(tr1))
        continue;
      registry.get<TH1>(HIST("Reco/trpt_k"))->Fill(vr0.Pt(), 1.);
      int t2 = 0;
      for (const auto& tr2 : tracks) {
        if (!tr2.has_udMcParticle())
          continue;
        if (!trackselector(tr2, parameters))
          continue;
        t2++;
        if (t2 > t1) {
          if (!selectionPIDKaon1(tr2))
            continue;
          vr1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
          auto mcPart2 = tr2.udMcParticle();
          //  if (std::abs(mcPart2.globalIndex() - mcPart1.globalIndex()) != 1)
          //   continue;
          if (std::abs(mcPart1.pdgCode()) != kKPlus || std::abs(mcPart2.pdgCode()) != kKPlus)
            continue;
          if (mcPart1.pdgCode() == mcPart2.pdgCode())
            continue;
          if (!mcPart1.isPhysicalPrimary() || !mcPart2.isPhysicalPrimary())
            continue;
          bool flag = false;
          bool flag1 = false;
          int gIndex1 = 0;
          int gIndex2 = 0;
          if (mcPart1.has_mothers() && mcPart2.has_mothers()) {
            for (const auto& mother : mcPart1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
                if (mother.has_mothers()) {
                  continue;
                }
                vphi.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassPhi);
                if (std::abs(vphi.Rapidity()) > 0.5)
                  continue;
                flag = true;
                gIndex1 = mother.globalIndex();
              }
            }
            for (const auto& mother1 : mcPart1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
                vphi.SetCoordinates(mother1.px(), mother1.py(), mother1.pz(), o2::constants::physics::MassPhi);
                if (std::abs(vphi.Rapidity()) > 0.5)
                  continue;
                flag1 = true;
                gIndex2 = mother1.globalIndex();
              }
            }
          }
          vr0g.SetCoordinates(mcPart1.px(), mcPart1.py(), mcPart1.pz(), o2::constants::physics::MassKaonCharged);
          vr1g.SetCoordinates(mcPart2.px(), mcPart2.py(), mcPart2.pz(), o2::constants::physics::MassKaonCharged);
          vr01g = vr0g + vr1g;
          vr01 = vr0 + vr1;
          if (flag && flag1 && (gIndex1 == gIndex2)) {
            if (std::abs(vr01g.Rapidity()) < 0.5) {
              registry.get<TH2>(HIST("Reco/selMPt"))->Fill(vr01g.M(), vr01g.Pt(), 1.);
            }
            if (std::abs(vr01.Rapidity()) < 0.5) {
              registry.get<TH3>(HIST("Reco/selMPtRap"))->Fill(vr01.M(), vr01.Pt(), vr01.Rapidity(), 1.);
            }
            if (std::abs(vphi.Rapidity()) < 0.5) {
              registry.get<TH3>(HIST("MC/accMPtRap_phi_T"))->Fill(vphi.M(), vphi.Pt(), vphi.Rapidity(), 1.);
            }
          }
          registry.get<TH1>(HIST("Reco/selRap"))->Fill(vr01.Rapidity(), 1.);
          registry.get<TH1>(HIST("Reco/selPt"))->Fill(vr01.Pt(), 1.);
          registry.get<TH1>(HIST("Reco/selM"))->Fill(vr01.M(), 1.);
          // registry.get<TH3>(HIST("Reco/selMPtRap_gen"))->Fill(vr01g.M(), vr01g.Pt(), vr01g.Rapidity(), 1.);
        }
      }
    }
    // KStar
    for (const auto& [tr1, tr2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!tr1.has_udMcParticle() || !tr2.has_udMcParticle())
        continue;
      if (!selectionPIDPion1(tr1) || !selectionPIDKaon1(tr2))
        continue;
      //  if (tr1.index() == tr2.index())
      // continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.
      auto mcPart1 = tr1.udMcParticle();
      auto mcPart2 = tr2.udMcParticle();
      if (std::abs(mcPart1.pdgCode()) != kPiPlus || std::abs(mcPart2.pdgCode()) != kKPlus)
        continue;
      if (!mcPart1.isPhysicalPrimary() || !mcPart2.isPhysicalPrimary())
        continue;
      // if (std::abs(mcPart2.globalIndex() - mcPart1.globalIndex()) != 1)
      //  continue;

      if (tr1.sign() * tr2.sign() > 0)
        continue;

      vr0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassPionCharged);
      vr1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
      vr0g.SetCoordinates(mcPart1.px(), mcPart1.py(), mcPart1.pz(), o2::constants::physics::MassPionCharged);
      vr1g.SetCoordinates(mcPart2.px(), mcPart2.py(), mcPart2.pz(), o2::constants::physics::MassKaonCharged);
      vr01g = vr0g + vr1g;
      vr01 = vr0 + vr1;
      if (!trackselector(tr1, parameters) || !trackselector(tr2, parameters)) {
        registry.get<TH1>(HIST("Reco/selM_k"))->Fill(vr01.M(), 1.);
      }
      if (trackselector(tr1, parameters) && trackselector(tr2, parameters)) {
        bool flag = false;
        bool flag1 = false;
        int gIndex1 = 0;
        int gIndex2 = 0;
        if (mcPart1.has_mothers() && mcPart2.has_mothers()) {
          for (const auto& mother : mcPart1.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
              vkstar.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassK0Star892);
              // registry.get<TH3>(HIST("MC/accMPtRap_kstar_T"))->Fill(vkstar.M(), vkstar.Pt(), vkstar.Rapidity(), 1.);
              if (std::abs(vkstar.Rapidity()) > 0.5)
                continue;
              flag = true;
              gIndex1 = mother.globalIndex();
            }
          }
          for (const auto& mother1 : mcPart2.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
              vkstar.SetCoordinates(mother1.px(), mother1.py(), mother1.pz(), o2::constants::physics::MassK0Star892);
              if (std::abs(vkstar.Rapidity()) > 0.5)
                continue;
              flag1 = true;
              gIndex2 = mother1.globalIndex();
            }
          }
        }
        if (flag && flag1 && (gIndex1 == gIndex2)) {
          if (std::abs(vr01g.Rapidity()) < 0.5) {
            registry.get<TH2>(HIST("Reco/selMPt_k"))->Fill(vr01g.M(), vr01g.Pt(), 1.);
          }
          if (std::abs(vr01.Rapidity()) < 0.5) {
            registry.get<TH3>(HIST("Reco/selMPtRap_k"))->Fill(vr01.M(), vr01.Pt(), vr01.Rapidity(), 1.);
          }
          if (std::abs(vkstar.Rapidity()) < 0.5) {
            registry.get<TH3>(HIST("MC/accMPtRap_kstar_T"))->Fill(vkstar.M(), vkstar.Pt(), vkstar.Rapidity(), 1.);
          }
        }
        registry.get<TH1>(HIST("Reco/selM_k_K"))->Fill(vr01.M(), 1.);
        registry.get<TH1>(HIST("Reco/selRap_k"))->Fill(vr01.Rapidity(), 1.);

        registry.get<TH1>(HIST("Reco/selPt_k"))->Fill(vr01.Pt(), 1.);
        registry.get<TH3>(HIST("Reco/selMPtRap_k_gen"))->Fill(vr01g.M(), vr01g.Pt(), vr01g.Rapidity(), 1.);
      }
    }
    registry.get<TH1>(HIST("Reco/nTracks"))->Fill(t1, 1.);
    // now access the McTruth information
    // get McCollision belonging to collision
    if (collision.has_udMcCollision()) {
      // auto mccollision = collision.udMcCollision();
      registry.get<TH1>(HIST("Reco/Stat"))->Fill(3., 1.);
    }
    // compute the difference between generated and reconstructed momentum
    for (const auto& track : tracks) {
      // is there an associated McParticle?
      if (track.has_udMcParticle()) {
        auto pTrack = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
        auto mcPart = track.udMcParticle();
        auto pPart = std::sqrt(mcPart.px() * mcPart.px() + mcPart.py() * mcPart.py() + mcPart.pz() * mcPart.pz());
        auto pDiff = pTrack - pPart;
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
      } else {
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(-5.9, -1, 1.);
      }
    }
  }
  PROCESS_SWITCH(SginclusivePhiKstarSD, processReco, "Process reconstructed data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SginclusivePhiKstarSD>(cfgc)};
}
