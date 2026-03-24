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

///
/// \file   qaEventTrackLite.cxx
/// \author Mario Krüger <mario.kruger@cern.ch>
/// \author Mattia Faggin <mattia.faggin@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \brief  Light version of the task to produce QA objects for the track and the event properties.
///         This task runs on prefiltered data
///

#include "qaEventTrack.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "MathUtils/BetheBlochAleph.h"
#include "ReconstructionDataFormats/PID.h"

#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

// Lite version of the QA task to run on skimmed dataset
struct qaEventTrackLite {
  // Binning
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis bins1overPt{"bins1overPt", {100, 0, 10}, "1/pt binning (c/GeV)"};
  ConfigurableAxis binsImpPar{"binsImpPar", {200, -0.15, 0.15}, "Impact parameter binning (cm)"};
  ConfigurableAxis binsEta{"binsEta", {200, -2., 2.}, "Eta binning"};
  ConfigurableAxis binsPhi{"binsPhi", {180, 0., 2 * M_PI}, "Phi binning"};
  ConfigurableAxis binsVtxZ{"binsVtxZ", {100, -20., 20.}, "Vtx Z binning"};

  HistogramRegistry histos;

  // Event selections
  Configurable<float> vtxZMax{"vtxZMax", 10.f, "Max VTX Z position"};
  // Track selections
  Configurable<bool> bHasITS{"bHasITS", false, "Select only DPG tracks with ITS (no forcing on TPC information)"};
  Configurable<bool> bItsStandalone{"bItsStandalone", false, "Select only ITS standalone DPG tracks"};
  Configurable<bool> bTpcOnly{"bTpcOnly", false, "Select only TPC only DPG tracks"};
  Configurable<bool> bItsTpcMatched{"bItsTpcMatched", false, "Select ITS-TPC matched DPG tracks"};
  // Kinematic selections
  Configurable<float> ptMin{"ptMin", 0., "Minimum track pt"};
  Configurable<float> etaMin{"etaMin", -10., "Minimum eta for DPG tracks"};
  Configurable<float> etaMax{"etaMax", 10., "Maximum eta for DPG tracks"};
  // ITS selections
  Configurable<float> chi2ItsMin{"chi2ItsMin", -1001.f, "Max ITS chi2"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 1000.f, "Max ITS chi2"};
  // TPC selections
  Configurable<int> nClusterTpcMin{"nClusterTpcMin", -1001, "Minimum number of TPC clusters"};
  Configurable<int> nCrossedRowsTpcMin{"nCrossedRowsTpcMin", -1001, "Minimum number of TPC crossed rows"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMin{"nCrossedRowsTpcOverFindableClustersTpcMin", -1, "Minimum ratio between TPC crossed rows and findable clusters"};
  Configurable<float> chi2TpcMin{"chi2TpcMin", -1001.f, "Max TPC chi2"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 1000.f, "Max TPC chi2"};
  Configurable<std::string> tpcSplinesPeriod{"tpcSplinesPeriod", std::string(""), "Period of used TPC dEdx splines"};
  Configurable<bool> b_tpcResProton{"b_tpcResProton", false, "Do TPC dEdx residuals around proton hypothesis"};
  Configurable<bool> b_tpcResKaon{"b_tpcResKaon", false, "Do TPC dEdx residuals around kaon hypothesis"};
  Configurable<bool> b_tpcResPion{"b_tpcResPion", false, "Do TPC dEdx residuals around pion hypothesis"};
  // TOF selections
  Configurable<float> chi2TofMin{"chi2TofMin", -1001.f, "Max TOF chi2"};
  Configurable<float> lengthMin{"lengthMin", -1001.f, "Min length"};
  Configurable<float> dcaXYmax{"dcaXYMax", 999., "Max dca XY"};

  // Selection 1
  // ITS selections
  Configurable<float> chi2ItsMaxSel1{"chi2ItsMaxSel1", 1000.f, "Max ITS chi2 Sel1"};
  // TPC selections
  Configurable<int> nClusterTpcMinSel1{"nClusterTpcMinSel1", -1001, "Minimum number of TPC clusters Sel1"};
  Configurable<int> nCrossedRowsTpcMinSel1{"nCrossedRowsTpcMinSel1", -1001, "Minimum number of TPC crossed rows Sel1"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMinSel1{"nCrossedRowsTpcOverFindableClustersTpcMinSel1", -1, "Minimum ratio between TPC crossed rows and findable clusters Sel1"};
  Configurable<float> chi2TpcMaxSel1{"chi2TpcMaxSel1", 1000.f, "Max TPC chi2 Sel1"};
  Configurable<bool> bItsTpcMatchedSel1{"bItsTpcMatchedSel1", false, "Select ITS-TPC matched DPG tracks, sel1"};
  Configurable<float> dcaXYmaxSel1{"dcaXYMaxSel1", 999., "Max dca XY sel1"};

  // Selection 2
  // ITS selections
  Configurable<float> chi2ItsMaxSel2{"chi2ItsMaxSel2", 1000.f, "Max ITS chi2 Sel2"};
  // TPC selections
  Configurable<int> nClusterTpcMinSel2{"nClusterTpcMinSel2", -1001, "Minimum number of TPC clusters Sel2"};
  Configurable<int> nCrossedRowsTpcMinSel2{"nCrossedRowsTpcMinSel2", -1001, "Minimum number of TPC crossed rows Sel2"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMinSel2{"nCrossedRowsTpcOverFindableClustersTpcMinSel2", -1, "Minimum ratio between TPC crossed rows and findable clusters Sel2"};
  Configurable<float> chi2TpcMaxSel2{"chi2TpcMaxSel2", 1000.f, "Max TPC chi2 Sel2"};
  Configurable<bool> bItsTpcMatchedSel2{"bItsTpcMatchedSel2", false, "Select ITS-TPC matched DPG tracks, sel2"};
  Configurable<float> dcaXYmaxSel2{"dcaXYMaxSel2", 999., "Max dca XY sel2"};

  // Selection 3
  // ITS selections
  Configurable<float> chi2ItsMaxSel3{"chi2ItsMaxSel3", 1000.f, "Max ITS chi2 Sel3"};
  // TPC selections
  Configurable<int> nClusterTpcMinSel3{"nClusterTpcMinSel3", -1001, "Minimum number of TPC clusters Sel3"};
  Configurable<int> nCrossedRowsTpcMinSel3{"nCrossedRowsTpcMinSel3", -1001, "Minimum number of TPC crossed rows Sel3"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMinSel3{"nCrossedRowsTpcOverFindableClustersTpcMinSel3", -1, "Minimum ratio between TPC crossed rows and findable clusters Sel3"};
  Configurable<float> chi2TpcMaxSel3{"chi2TpcMaxSel3", 1000.f, "Max TPC chi2 Sel3"};
  Configurable<bool> bItsTpcMatchedSel3{"bItsTpcMatchedSel3", false, "Select ITS-TPC matched DPG tracks, sel3"};
  Configurable<float> dcaXYmaxSel3{"dcaXYMaxSel3", 999., "Max dca XY sel3"};

  // MC selections
  Configurable<int> pdgCodeSel{"pdgCodeSel", 0, "pdgCode based particle selection. Provide a PDG code required for particles to have. To be used in combo with pdgCodeSelMode"};
  Configurable<int> pdgCodeSelMode{"pdgCodeSelMode", 2, "multiple pdgCode based particle selection. `1` accepts pi,K,p,mu,e, `2` accepts all final-state charged particles including light (hyper)nuclei"};

  Configurable<bool> checkPdgAtReco{"checkPdgAtReco", false, "check pdg code also at reco levo for data-like reference"};

  // TPC dEdx splines
  struct BBAleph {

    // --- data members ---
    double mMip = 0;
    std::vector<double> mBetheBlockAleph = {0., 0., 0., 0., 0.};
    double mChargeFactor = 0;
    bool initBBok = false;

    // --- functions ---
    double operator()(double* x, double* par)
    {
      /// === Parameters ===
      ///  [0]: mass
      ///  [1]: charge
      ///
      /// From A. Kalteyer:
      /// const float bethe = mMIP
      ///                   * o2::common::BetheBlochAleph(track.tpcInnerParam() / o2::track::pid_constants::sMasses[id]
      ///                   , mBetheBlochParams[0], mBetheBlochParams[1], mBetheBlochParams[2], mBetheBlochParams[3], mBetheBlochParams[4])
      ///                   * std::pow((float)o2::track::pid_constants::sCharges[id], mChargeFactor);
      ///
      return initBBok ? mMip * o2::common::BetheBlochAleph(x[0] / par[0], mBetheBlockAleph[0], mBetheBlockAleph[1], mBetheBlockAleph[2], mBetheBlockAleph[3], mBetheBlockAleph[4]) * std::pow(par[1], mChargeFactor) : 0.;
    }
    void setUpBetheBlockAleph(std::string str_case)
    {
      if (str_case.find("LHC22c") != std::string::npos) {
        // From A. Kalteyer (2022 Jul 18)
        mMip = 52.35295;
        mChargeFactor = 4.382516;
        mBetheBlockAleph[0] = 0.038508328324810444;
        mBetheBlockAleph[1] = 20.173734667349745;
        mBetheBlockAleph[2] = 7.811839314098901e-10;
        mBetheBlockAleph[3] = 2.3572347501513162;
        mBetheBlockAleph[4] = 3.896773574265881;
        initBBok = true;
      } else if (str_case.find("LHC22d") != std::string::npos) {
        // From A. Kalteyer (2022 Jul 18)
        mMip = 52.35295;
        mChargeFactor = 4.382516;
        mBetheBlockAleph[0] = 0.0370972540453123;
        mBetheBlockAleph[1] = 21.29155104648775;
        mBetheBlockAleph[2] = 1.2702695962030138e-10;
        mBetheBlockAleph[3] = 2.3297081952872936;
        mBetheBlockAleph[4] = 4.1832594396288725;
        initBBok = true;
      } else {
        LOG(info) << "===> WARNING: no Bethe Block parameters defined for " << str_case << ". Ignoring it." << std::endl;
        return;
      }
      LOG(info) << "=== Bethe Block Alep parametrization loaded for period " << str_case << ":" << std::endl;
      LOG(info) << "mMip                = " << mMip << std::endl;
      LOG(info) << "mChargeFactor       = " << mChargeFactor << std::endl;
      LOG(info) << "mBetheBlockAleph[0] = " << mBetheBlockAleph[0] << std::endl;
      LOG(info) << "mBetheBlockAleph[1] = " << mBetheBlockAleph[1] << std::endl;
      LOG(info) << "mBetheBlockAleph[2] = " << mBetheBlockAleph[2] << std::endl;
      LOG(info) << "mBetheBlockAleph[3] = " << mBetheBlockAleph[3] << std::endl;
      LOG(info) << "mBetheBlockAleph[4] = " << mBetheBlockAleph[4] << std::endl;
      LOG(info) << "initBBok            = " << initBBok << std::endl;
    }
  };
  BBAleph betheBlock;
  TF1 funcBBpion, funcBBkaon, funcBBproton; // TPC dEdx splines

  void init(InitContext const&)
  {
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axis1overPt{bins1overPt, "1/#it{p}_{T} (GeV/c)^{-1}"};
    const AxisSpec axisEta{binsEta, "#it{#eta}"};
    const AxisSpec axisPhi{binsPhi, "#it{#phi} (rad)"};
    const AxisSpec axisVtxZ{binsVtxZ, "Vertex Z (cm)"};

    // TPC dEdx splines
    betheBlock.setUpBetheBlockAleph(tpcSplinesPeriod);

    // kine histograms
    histos.add("Tracks/VertexPositionZ", "", kTH1D, {axisVtxZ});
    histos.add("Tracks/VertexPositionZvsEta", "", kTH2D, {axisEta, axisVtxZ});
    histos.add("Tracks/VertexPositionZvsEtaHasITS", "", kTH2D, {axisEta, axisVtxZ});
    histos.add("Tracks/Kine/pt", "#it{p}_{T}", kTH1D, {axisPt});
    histos.add("Tracks/Kine/eta", "#eta", kTH1D, {axisEta});
    histos.add("Tracks/Kine/phi", "#phi", kTH1D, {axisPhi});
    histos.add("Tracks/Kine/phiVsEta", "#phi", kTH2D, {axisPhi, axisEta});
    histos.add("Tracks/length", "track length in cm;#it{Length} (cm);", kTH1D, {{400, 0, 1000}});
    const AxisSpec axisImpParRPhi{binsImpPar, "#it{d}_{r#it{#varphi}} (#cm)"};
    const AxisSpec axisImpParZAxis{binsImpPar, "#it{d}_{z} (#cm)"};
    histos.add("Tracks/dcaXY", "distance of closest approach in #it{xy} plane", kTH1D, {axisImpParRPhi});
    histos.add("Tracks/dcaZ", "distance of closest approach in #it{z}", kTH1D, {axisImpParZAxis});
    histos.add("Tracks/dcaXYvsPt", "d_#it{xy} vs. #it{p}_{T}", kTH2D, {axisPt, axisImpParRPhi});
    histos.add("Tracks/dcaZvsPt", "d_#it{z} vs. #it{p}_{T}", kTH2D, {axisPt, axisImpParRPhi});
    histos.add("Tracks/dcaXYvsPtvsEta", "d_#it{xy} vs. #it{p}_{T} vs. #eta", kTH3D, {axisPt, axisEta, axisImpParRPhi});
    histos.add("Tracks/dcaZvsPtvsEta", "d_#it{z} vs. #it{p}_{T} vs. #eta", kTH3D, {axisPt, axisEta, axisImpParRPhi});

    // its histograms
    histos.add("Tracks/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {{100, 0, 40}});
    histos.add("Tracks/ITS/itsNCl", "ITS number of clusters;# clusters ITS", kTH1D, {{8, -0.5, 7.5}});
    histos.add("Tracks/ITS/itsNClvsItsHitmap", "ITS number of clusters vs. ITS hitmap;# clusters ITS; ITS hitmap", kTH2D, {{8, -0.5, 7.5, "# clusters ITS"}, {128, 0, 128, "ITS hitmap"}});
    histos.add("Tracks/ITS/itsNClstvsEtavsPt", "profile2D;", kTProfile2D, {axisEta, axisPt});
    // tpc histograms
    histos.add("Tracks/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {{100, 0, 10}});
    histos.add("Tracks/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {{165, -0.5, 164.5}});
    histos.add("Tracks/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {{165, -0.5, 164.5}});
    histos.add("Tracks/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {{60, 0.7, 1.3}});
    histos.add("Tracks/TPC/tpcNClsFoundvsPt", "", kTH2D, {axisPt, {165, -0.5, 164.5, "# clusters TPC"}});
    histos.add("Tracks/TPC/tpcCrossedRowsvsPt", "", kTH2D, {axisPt, {165, -0.5, 164.5, "# crossed rows TPC"}});
    histos.add("Tracks/TPC/tpcCrossedRowsOverFindableClsvsPt", "", kTH2D, {axisPt, {60, 0.7, 1.3, "crossed rows / findable clusters TPC"}});
    histos.add("Tracks/TPC/TPCnClstvsEtavsPt", "profile2D;", kTProfile2D, {axisEta, axisPt});
    histos.add("Tracks/TPC/dEdxvsP", "", kTH2D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {500, 0, 1000, "d#it{E}/d#it{x} (a.u.)"}});
    histos.add("Tracks/TPC/dEdxvsPvsEta", "", kTH3D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {20, -2, 2, "#it{#eta}"}, {500, 0, 1000, "d#it{E}/d#it{x} (a.u.)"}});
    if (betheBlock.initBBok) {
      if (b_tpcResProton) {
        // Bethe-Block parametrization proton
        funcBBproton = TF1("funcBBproton", betheBlock, 0.1, 20, 3);
        funcBBproton.SetParameter(0, o2::track::pid_constants::sMasses[o2::track::PID::Proton]);  // mass
        funcBBproton.SetParameter(1, o2::track::pid_constants::sCharges[o2::track::PID::Proton]); // electric charge (units of e)
        // histograms of residuals
        histos.add("Tracks/TPC/dEdxvsPproton", "", kTH2D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {500, -100, 100, "d#it{E}/d#it{x}-d#it{E}/d#it{x}|_{proton} (a.u.)"}});
        histos.add("Tracks/TPC/dEdxvsPprotonvsEta", "", kTH3D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {20, -2, 2, "#it{#eta}"}, {500, -100, 100, "d#it{E}/d#it{x}-d#it{E}/d#it{x}|_{proton} (a.u.)"}});
      }
      if (b_tpcResKaon) {
        // Bethe-Block parametrization kaon
        funcBBkaon = TF1("funcBBkaon", betheBlock, 0.1, 20, 3);
        funcBBkaon.SetParameter(0, o2::track::pid_constants::sMasses[o2::track::PID::Kaon]);  // mass
        funcBBkaon.SetParameter(1, o2::track::pid_constants::sCharges[o2::track::PID::Kaon]); // electric charge (units of e)
        // histograms of residuals
        histos.add("Tracks/TPC/dEdxvsPkaon", "", kTH2D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {500, -100, 100, "d#it{E}/d#it{x}-d#it{E}/d#it{x}|_{kaon} (a.u.)"}});
        histos.add("Tracks/TPC/dEdxvsPkaonvsEta", "", kTH3D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {20, -2, 2, "#it{#eta}"}, {500, -100, 100, "d#it{E}/d#it{x}-d#it{E}/d#it{x}|_{kaon} (a.u.)"}});
      }
      if (b_tpcResPion) {
        // Bethe-Block parametrization pion
        funcBBpion = TF1("funcBBpion", betheBlock, 0.1, 20, 3);
        funcBBpion.SetParameter(0, o2::track::pid_constants::sMasses[o2::track::PID::Pion]);  // mass
        funcBBpion.SetParameter(1, o2::track::pid_constants::sCharges[o2::track::PID::Pion]); // electric charge (units of e)
        // histograms of residuals
        histos.add("Tracks/TPC/dEdxvsPpion", "", kTH2D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {500, -100, 100, "d#it{E}/d#it{x}-d#it{E}/d#it{x}|_{pion} (a.u.)"}});
        histos.add("Tracks/TPC/dEdxvsPpionvsEta", "", kTH3D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {20, -2, 2, "#it{#eta}"}, {500, -100, 100, "d#it{E}/d#it{x}-d#it{E}/d#it{x}|_{pion} (a.u.)"}});
      }
    }
    // trd histograms
    histos.add("Tracks/TRD/trdChi2", "chi2 in TRD", kTH1D, {{100, 0, 10, "chi2 / cluster TRD"}});
    // tof histograms
    histos.add("Tracks/TOF/tofChi2", "chi2 in TOF", kTH1D, {{100, 0, 10, "chi2 / cluster TOF"}});
    // matching histogram
    histos.add("Tracks/matchedDet", "matched detectors", kTH1D, {{4, 0.5, 4.5, ""}});
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(1, "hasTPC");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(2, "hasITS");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(3, "hasTRD");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(4, "hasTOF");
    // kinematics
    histos.add("Tracks/relativeResoPt", "relative #it{p}_{T} resolution;#sigma(#it{p}_{T})/#it{p}_{T};#it{p}_{T};", kTH2D, {{axisPt, {500, 0., 1, "#sigma{#it{p}}/#it{p}_{T}"}}});
    histos.add("Tracks/relativeResoPtMean", "mean relative #it{p}_{T} resolution;;average #sigma_{#it{p}_{T}}/#it{p}_{T};", kTProfile, {axisPt});
    histos.add("Tracks/relativeResoPtMeanvsEtavsPt", "mean relative #it{p}_{T} resolution;", kTProfile2D, {axisPt, axisEta});
    histos.add("Tracks/reso1overPtMeanvsEtavs1overPt", "mean 1/#it{p}_{T} resolution;", kTProfile2D, {axis1overPt, axisEta});
    // track cuts map
    histos.add("TrackCuts/selPtEtaPhiNoSel", "pt eta phi map no sel; pt,eta,phi", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("TrackCuts/selPtEtaPhiSel1", "pt eta phi map sel 1; pt,eta,phi", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("TrackCuts/selPtEtaPhiSel2", "pt eta phi map sel 2; pt,eta,phi", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("TrackCuts/selPtEtaPhiSel3", "pt eta phi map sel 3; pt,eta,phi", kTH3D, {axisPt, axisEta, axisPhi});

    // MC histograms
    if (doprocessMCLite) {
      histos.add("Particles/PDGs", "Particle PDGs;PDG Code", kTH1D, {{100, 0.f, 100.f}});
      histos.add("Particles/Kine/pt", "Particle #it{p}_{T}", kTH1D, {axisPt});
      histos.add("Particles/Kine/eta", "Particle #eta", kTH1D, {axisEta});
      histos.add("Particles/Kine/phi", "Particle #phi", kTH1D, {axisPhi});
      histos.add("Particle/selPtEtaPhiMCGenPrimary", "pt eta phi map MC gen Primary", kTH3D, {axisPt, axisEta, axisPhi});
      histos.add("Particle/selPtEtaPhiMCRecoNoSelPrimary", "pt eta phi map MC reco Primary sel0", kTH3D, {axisPt, axisEta, axisPhi});
      histos.add("Particle/selPtEtaPhiMCRecoSel1Primary", "pt eta phi map MC reco Primary sel1", kTH3D, {axisPt, axisEta, axisPhi});
      histos.add("Particle/selPtEtaPhiMCRecoSel2Primary", "pt eta phi map MC reco Primary sel2", kTH3D, {axisPt, axisEta, axisPhi});
      histos.add("Particle/selPtEtaPhiMCRecoSel3Primary", "pt eta phi map MC reco Primary sel3", kTH3D, {axisPt, axisEta, axisPhi});

      histos.add("Particle/selPtEtaVtxzMCGenPrimary", "pt eta VtxZ map MC gen Primary", kTH3F, {axisPt, axisEta, axisVtxZ});
      histos.add("Particle/selPtEtaVtxzMCRecoNoSelPrimary", "pt eta VtxZ map MC reco Primary sel0", kTH3F, {axisPt, axisEta, axisVtxZ});
      histos.add("Particle/selPtEtaVtxzMCRecoSel1Primary", "pt eta VtxZ map MC reco Primary sel1", kTH3F, {axisPt, axisEta, axisVtxZ});
      histos.add("Particle/selPtEtaVtxzMCRecoSel2Primary", "pt eta VtxZ map MC reco Primary sel2", kTH3F, {axisPt, axisEta, axisVtxZ});
      histos.add("Particle/selPtEtaVtxzMCRecoSel3Primary", "pt eta VtxZ map MC reco Primary sel3", kTH3F, {axisPt, axisEta, axisVtxZ});
      histos.add("Tracks/MC/resoPhivsPtvsEta", "#varphi(reco)-#varphi(gen);", kTH3D, {axisPt, axisEta, {180, -M_PI, M_PI, "#varphi(reco)-#varphi(gen)"}});
      histos.add("Tracks/MC/phiRecovsphiGen", "#varphi(reco) vs. #varphi(gen);", kTH2D, {axisPhi, axisPhi});
      histos.get<TH2>(HIST("Tracks/MC/phiRecovsphiGen"))->GetXaxis()->SetTitle("#varphi(reco)");
      histos.get<TH2>(HIST("Tracks/MC/phiRecovsphiGen"))->GetYaxis()->SetTitle("#varphi(gen)");
      histos.add("Tracks/MC/ptRecoVsptGen", "", kTH2D, {{axisPt, axisPt}});
      histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen"))->GetYaxis()->SetTitle(Form("%s_{Gen}", histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen"))->GetYaxis()->GetTitle()));
      histos.add("Tracks/MC/ptRecoVsptGen_wTOF", "", kTH2D, {{axisPt, axisPt}});
      histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen_wTOF"))->GetYaxis()->SetTitle(Form("%s_{Gen}", histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen_wTOF"))->GetYaxis()->GetTitle()));
      histos.add("Tracks/MC/ptRecoVsptGen_wTRD", "", kTH2D, {{axisPt, axisPt}});
      histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen_wTRD"))->GetYaxis()->SetTitle(Form("%s_{Gen}", histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen_wTRD"))->GetYaxis()->GetTitle()));
      histos.add("Tracks/MC/ptRecoVsptGen_woTRD", "", kTH2D, {{axisPt, axisPt}});
      histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen_woTRD"))->GetYaxis()->SetTitle(Form("%s_{Gen}", histos.get<TH2>(HIST("Tracks/MC/ptRecoVsptGen_woTRD"))->GetYaxis()->GetTitle()));
    }
  }

  ///////////////
  /// Filters ///
  ///////////////
  // Kinematics
  Filter ptCut = o2::aod::dpgtrack::pt > ptMin;
  Filter etaCut = etaMin < o2::aod::dpgtrack::eta && o2::aod::dpgtrack::eta < etaMax;
  // Detector matching
  Filter filterHasITS = (bHasITS.node() == false) || (bItsStandalone.node() == false && bTpcOnly.node() == false && bItsTpcMatched.node() == false && o2::aod::dpgtrack::hasITS == true);
  Filter itsStandaloneTracks = (bItsStandalone.node() == false) || (o2::aod::dpgtrack::hasITS == true && o2::aod::dpgtrack::hasTPC == false);
  Filter tpcOnlyTracks = (bTpcOnly.node() == false) || (o2::aod::dpgtrack::hasITS == false && o2::aod::dpgtrack::hasTPC == true);
  Filter itsTpcMatchedTracks = (bItsTpcMatched.node() == false) || (o2::aod::dpgtrack::hasITS == true && o2::aod::dpgtrack::hasTPC == true);
  // ITS
  Filter itsChi2 = (bTpcOnly.node() == true) || (chi2ItsMin < o2::aod::track::itsChi2NCl && o2::aod::track::itsChi2NCl < chi2ItsMax);
  // TPC
  Filter tpcChi2s = (bItsStandalone.node() == true) || (chi2TpcMin < o2::aod::track::tpcChi2NCl && o2::aod::track::tpcChi2NCl < chi2TpcMax);
  Filter tpcNclusters = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcNClsFound > (int16_t)nClusterTpcMin);
  Filter tpcNcrossedRows = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcNClsCrossedRows > (int16_t)nCrossedRowsTpcMin);
  Filter tpcNcrossedRowsOverFindableClusters = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcCrossedRowsOverFindableCls > nCrossedRowsTpcOverFindableClustersTpcMin);
  // TOF
  Filter tofChi = o2::aod::track::tofChi2 > chi2TofMin;
  Filter length = o2::aod::track::length > lengthMin;
  template <bool isMC, typename trackType>
  void fillHistograms(trackType const& track)
  {
    const float vtxZ = track.dpgCollision().posZ();
    if (TMath::Abs(vtxZ) > vtxZMax) {
      return;
    }
    if constexpr (isMC) {
      if (track.isPhysicalPrimary() && isPdgSelected(track.pdgCode())) {
        histos.fill(HIST("Particle/selPtEtaPhiMCGenPrimary"), track.ptMC(), track.etaMC(), track.phiMC());
        histos.fill(HIST("Particle/selPtEtaVtxzMCGenPrimary"), track.ptMC(), track.etaMC(), vtxZ);
      }
    }
    // temporary additional selections
    if (!applyTrackSelectionsNotFiltered(track)) {
      return;
    }
    Bool_t sel1 = applyTrackSelectionsSel1(track);
    Bool_t sel2 = applyTrackSelectionsSel2(track);
    Bool_t sel3 = applyTrackSelectionsSel3(track);
    Bool_t isPdgOk = true;
    if constexpr (isMC) {
      isPdgOk = isPdgSelected(track.pdgCode());
      if (track.isPhysicalPrimary() && isPdgOk) {
        histos.get<TH1>(HIST("Particles/PDGs"))->Fill(Form("%i", track.pdgCode()), 1);
        histos.fill(HIST("Particle/selPtEtaPhiMCRecoNoSelPrimary"), track.ptMC(), track.etaMC(), track.phiMC());
        histos.fill(HIST("Particle/selPtEtaVtxzMCRecoNoSelPrimary"), track.ptMC(), track.etaMC(), vtxZ);

        if (sel1) {
          histos.fill(HIST("Particle/selPtEtaPhiMCRecoSel1Primary"), track.ptMC(), track.etaMC(), track.phiMC());
          histos.fill(HIST("Particle/selPtEtaVtxzMCRecoSel1Primary"), track.ptMC(), track.etaMC(), vtxZ);
        }
        if (sel2) {
          histos.fill(HIST("Particle/selPtEtaPhiMCRecoSel2Primary"), track.ptMC(), track.etaMC(), track.phiMC());
          histos.fill(HIST("Particle/selPtEtaVtxzMCRecoSel2Primary"), track.ptMC(), track.etaMC(), vtxZ);
        }
        if (sel3) {
          histos.fill(HIST("Particle/selPtEtaPhiMCRecoSel3Primary"), track.ptMC(), track.etaMC(), track.phiMC());
          histos.fill(HIST("Particle/selPtEtaVtxzMCRecoSel3Primary"), track.ptMC(), track.etaMC(), vtxZ);
        }
      }

      histos.fill(HIST("Particles/Kine/pt"), track.ptMC());
      histos.fill(HIST("Particles/Kine/eta"), track.etaMC());
      histos.fill(HIST("Particles/Kine/phi"), track.phiMC());
    }
    histos.fill(HIST("Tracks/VertexPositionZ"), track.dpgCollision().posZ());
    histos.fill(HIST("Tracks/VertexPositionZvsEta"), track.eta(), track.dpgCollision().posZ());
    if (track.hasITS()) {
      histos.fill(HIST("Tracks/VertexPositionZvsEtaHasITS"), track.eta(), track.dpgCollision().posZ());
    }
    histos.fill(HIST("Tracks/Kine/pt"), track.pt());
    histos.fill(HIST("Tracks/Kine/eta"), track.eta());
    histos.fill(HIST("Tracks/Kine/phi"), track.phi());
    histos.fill(HIST("Tracks/Kine/phiVsEta"), track.phi(), track.eta());
    histos.fill(HIST("Tracks/dcaXY"), track.dcaXY());
    histos.fill(HIST("Tracks/dcaZ"), track.dcaZ());
    histos.fill(HIST("Tracks/dcaXYvsPt"), track.pt(), track.dcaXY());
    histos.fill(HIST("Tracks/dcaZvsPt"), track.pt(), track.dcaZ());
    histos.fill(HIST("Tracks/dcaXYvsPtvsEta"), track.pt(), track.eta(), track.dcaXY());
    histos.fill(HIST("Tracks/dcaZvsPtvsEta"), track.pt(), track.eta(), track.dcaZ());
    histos.fill(HIST("Tracks/length"), track.length());
    histos.fill(HIST("Tracks/ITS/itsChi2NCl"), track.itsChi2NCl());
    histos.fill(HIST("Tracks/ITS/itsNCl"), track.itsNCls());
    histos.fill(HIST("Tracks/ITS/itsNClvsItsHitmap"), track.itsNCls(), track.itsClusterMap());
    histos.fill(HIST("Tracks/TPC/tpcChi2NCl"), track.tpcChi2NCl());
    histos.fill(HIST("Tracks/TPC/tpcNClsFound"), track.tpcNClsFound());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
    histos.fill(HIST("Tracks/TPC/tpcNClsFoundvsPt"), track.pt(), track.tpcNClsFound());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRowsvsPt"), track.pt(), track.tpcNClsCrossedRows());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableClsvsPt"), track.pt(), track.tpcCrossedRowsOverFindableCls());
    const double p = track.pt() / (sin(2 * atan2(1, exp(track.eta()))));
    histos.fill(HIST("Tracks/TPC/dEdxvsP"), p, track.tpcSignal());
    histos.fill(HIST("Tracks/TPC/dEdxvsPvsEta"), p, track.eta(), track.tpcSignal());
    if (betheBlock.initBBok) {
      auto tpcdEdxRes = [&](TF1 func) { return track.tpcSignal() - func.Eval(p); };
      if (b_tpcResProton) {
        histos.fill(HIST("Tracks/TPC/dEdxvsPproton"), p, tpcdEdxRes(funcBBproton));
        histos.fill(HIST("Tracks/TPC/dEdxvsPprotonvsEta"), p, track.eta(), tpcdEdxRes(funcBBproton));
      }
      if (b_tpcResKaon) {
        histos.fill(HIST("Tracks/TPC/dEdxvsPkaon"), p, tpcdEdxRes(funcBBkaon));
        histos.fill(HIST("Tracks/TPC/dEdxvsPkaonvsEta"), p, track.eta(), tpcdEdxRes(funcBBkaon));
      }
      if (b_tpcResPion) {
        histos.fill(HIST("Tracks/TPC/dEdxvsPpion"), p, tpcdEdxRes(funcBBpion));
        histos.fill(HIST("Tracks/TPC/dEdxvsPpionvsEta"), p, track.eta(), tpcdEdxRes(funcBBpion));
      }
    }
    histos.fill(HIST("Tracks/TRD/trdChi2"), track.trdChi2());
    histos.fill(HIST("Tracks/TOF/tofChi2"), track.tofChi2());
    if (track.hasTPC()) {
      histos.fill(HIST("Tracks/matchedDet"), 1);
    }
    if (track.hasITS()) {
      histos.fill(HIST("Tracks/matchedDet"), 2);
    }
    if (track.hasTRD()) {
      histos.fill(HIST("Tracks/matchedDet"), 3);
    }
    if (track.hasTOF()) {
      histos.fill(HIST("Tracks/matchedDet"), 4);
    }
    histos.fill(HIST("Tracks/relativeResoPt"), track.pt(), track.ptReso());
    histos.fill(HIST("Tracks/relativeResoPtMean"), track.pt(), track.ptReso());
    histos.fill(HIST("Tracks/relativeResoPtMeanvsEtavsPt"), track.pt(), track.eta(), track.ptReso());
    histos.fill(HIST("Tracks/reso1overPtMeanvsEtavs1overPt"), 1. / track.pt(), track.eta(), track.ptReso() / (track.pt())); // reso(1/pt) = reso(pt)/pt^2 and here ptReso==reso(pt)/pt
    histos.fill(HIST("Tracks/TPC/TPCnClstvsEtavsPt"), track.eta(), track.pt(), track.tpcNClsFound());
    histos.fill(HIST("Tracks/ITS/itsNClstvsEtavsPt"), track.eta(), track.pt(), track.itsNCls());
    if constexpr (isMC) {
      histos.fill(HIST("Tracks/MC/ptRecoVsptGen"), track.pt(), track.ptMC());
      if (track.hasTOF()) {
        histos.fill(HIST("Tracks/MC/ptRecoVsptGen_wTOF"), track.pt(), track.ptMC());
      }
      if (track.hasTRD()) {
        histos.fill(HIST("Tracks/MC/ptRecoVsptGen_wTRD"), track.pt(), track.ptMC());
      } else {
        histos.fill(HIST("Tracks/MC/ptRecoVsptGen_woTRD"), track.pt(), track.ptMC());
      }
      histos.fill(HIST("Tracks/MC/resoPhivsPtvsEta"), track.pt(), track.eta(), track.phi() - track.phiMC());
      histos.fill(HIST("Tracks/MC/phiRecovsphiGen"), track.phi(), track.phiMC());
    }

    if constexpr (isMC) {
      if (isPdgOk || !checkPdgAtReco) {
        histos.fill(HIST("TrackCuts/selPtEtaPhiNoSel"), track.ptMC(), track.etaMC(), track.phiMC());
      }
    } else {
      histos.fill(HIST("TrackCuts/selPtEtaPhiNoSel"), track.pt(), track.eta(), track.phi());
    }

    if (sel1) {
      if constexpr (isMC) {
        if (isPdgOk || !checkPdgAtReco) {
          histos.fill(HIST("TrackCuts/selPtEtaPhiSel1"), track.ptMC(), track.etaMC(), track.phiMC());
        }
      } else {
        histos.fill(HIST("TrackCuts/selPtEtaPhiSel1"), track.pt(), track.eta(), track.phi());
      }
    }
    if (sel2) {
      if constexpr (isMC) {
        if (isPdgOk || !checkPdgAtReco) {
          histos.fill(HIST("TrackCuts/selPtEtaPhiSel2"), track.ptMC(), track.etaMC(), track.phiMC());
        }
      } else {
        histos.fill(HIST("TrackCuts/selPtEtaPhiSel2"), track.pt(), track.eta(), track.phi());
      }
    }
    if (sel3) {
      if constexpr (isMC) {
        if (isPdgOk || !checkPdgAtReco) {
          histos.fill(HIST("TrackCuts/selPtEtaPhiSel3"), track.ptMC(), track.etaMC(), track.phiMC());
        }
      } else {
        histos.fill(HIST("TrackCuts/selPtEtaPhiSel3"), track.pt(), track.eta(), track.phi());
      }
    }
  }

  // Process data
  void processDataLite(o2::soa::Filtered<aod::DPGTracks> const& tracks, aod::DPGCollisions const&)
  {
    for (const auto& track : tracks) {
      fillHistograms<false>(track);
    }
  }
  PROCESS_SWITCH(qaEventTrackLite, processDataLite, "process data lite", true);

  // Process MC
  void processMCLite(o2::soa::Filtered<soa::Join<aod::DPGTracks, aod::DPGRecoParticles>> const& tracks, aod::DPGCollisions const&,
                     aod::DPGNonRecoParticles const& nonRecoParticles)
  {
    for (const auto& track : tracks) {
      fillHistograms<true>(track);
    }

    for (const auto& particle : nonRecoParticles) {
      const float vtxZ = particle.dpgCollision().posZ();
      if (TMath::Abs(vtxZ) > vtxZMax) {
        continue;
      }

      histos.fill(HIST("Particles/Kine/pt"), particle.ptMC());
      histos.fill(HIST("Particles/Kine/eta"), particle.etaMC());
      histos.fill(HIST("Particles/Kine/phi"), particle.phiMC());

      if (particle.isPhysicalPrimary() && isPdgSelected(particle.pdgCode())) {
        histos.fill(HIST("Particle/selPtEtaPhiMCGenPrimary"), particle.ptMC(), particle.etaMC(), particle.phiMC());
        histos.fill(HIST("Particle/selPtEtaVtxzMCGenPrimary"), particle.ptMC(), particle.etaMC(), vtxZ);
      }
    }
  }
  PROCESS_SWITCH(qaEventTrackLite, processMCLite, "process MC lite", false);

  template <typename T>
  bool applyTrackSelectionsNotFiltered(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMin)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMin)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMin)
      return false;
    if (bItsStandalone == true && (track.hasITS() == false || track.hasTPC() == true))
      return false;
    if (bTpcOnly == true && (track.hasITS() == true || track.hasTPC() == false))
      return false;
    if (bItsTpcMatched == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.itsChi2NCl() > chi2ItsMax)
      return false;
    if (track.tpcChi2NCl() > chi2TpcMax)
      return false;
    if (track.tpcNClsFound() < nClusterTpcMin)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmax)
      return false;
    return true;
  }

  template <typename T>
  bool applyTrackSelectionsSel1(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMinSel1)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMinSel1)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMinSel1)
      return false;
    if (bItsTpcMatchedSel1 == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.itsChi2NCl() > chi2ItsMaxSel1)
      return false;
    if (track.tpcChi2NCl() > chi2TpcMaxSel1)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmaxSel1)
      return false;

    return true;
  }
  template <typename T>
  bool applyTrackSelectionsSel2(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMinSel2)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMinSel2)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMinSel2)
      return false;
    if (track.itsChi2NCl() > chi2ItsMaxSel2)
      return false;
    if (bItsTpcMatchedSel2 == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.tpcChi2NCl() > chi2TpcMaxSel2)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmaxSel2)
      return false;

    return true;
  }
  template <typename T>
  bool applyTrackSelectionsSel3(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMinSel3)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMinSel3)
      return false;
    if (bItsTpcMatchedSel3 == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMinSel3)
      return false;
    if (track.itsChi2NCl() > chi2ItsMaxSel3)
      return false;
    if (track.tpcChi2NCl() > chi2TpcMaxSel3)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmaxSel3)
      return false;

    return true;
  }

  bool isPdgSelected(const int pdgcode)
  { // mimics selection of charged particles or id particles

    if (pdgcode == pdgCodeSel) { // Check that the pdg code is exactly what was asked
      return true;
    } else if (pdgCodeSel != 0) {
      return false;
    }
    const int abspdgcode = abs(pdgcode);
    if (pdgCodeSelMode == 1 || pdgCodeSelMode == 2) {
      switch (abspdgcode) {
        case 11:   // electron
        case 13:   // muon
        case 211:  // pion
        case 321:  // kaon
        case 2212: // proton
          return true;
      }

      if (pdgCodeSelMode == 2) {
        switch (abspdgcode) {
          case 3222:       // Σ+
          case 3112:       // Σ−
          case 3312:       // Ξ−
          case 3334:       // Ω−
          case 1000010020: // deuteron
          case 1000010030: // triton
          case 1000020030: // helium3
          case 1000020040: // helium4
          case 1010010030: // hyper triton
          case 1010020040: // hyper helium4
            return true;
        }
      }
    }
    return false;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaEventTrackLite>(cfgc)};
}
