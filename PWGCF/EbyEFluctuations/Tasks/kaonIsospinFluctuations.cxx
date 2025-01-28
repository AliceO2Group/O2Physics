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

/// \file kaonIsospinFluctuations.cxx
/// \brief Kaon Isospin fluctuations
///
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Sadhana Dash (sadhana@phy.iitb.ac.in)

#include <algorithm>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics; // for constants

struct KaonIsospinFluctuations {
  // Hisogram registry:
  HistogramRegistry recoV0s{"recoV0s", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoEvent{"recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoK0s{"recoK0s", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoTracks{"recoTracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoAnalysis{"recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Configurables
  // Event Selection
  Configurable<float> cutZvertex{"cutZvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> v0settingDcaPosToPV{"v0settingDcaPosToPV", 0.06, "DCA Pos to PV"};
  Configurable<float> v0settingDcaNegToPV{"v0settingDcaNegToPV", 0.06, "DCA Neg to PV"};
  Configurable<float> v0settingDcaV0Dau{"v0settingDcaV0Dau", 1, "DCA V0 Daughters"};
  Configurable<double> v0settingCosPA{"v0settingCosPA", 0.98, "V0 CosPA"};
  Configurable<float> v0settingRadius{"v0settingRadius", 0.5, "v0radius"};

  // Configurable K0s  //0.480-0.515 GeV/c2
  Configurable<double> mLowK0s{"mLowK0s", 0.48, "mLowK0s"};
  Configurable<double> mHighK0s{"mHighK0s", 0.515, "mHighK0s"};

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec axisK0sMass = {200, 0.40f, 0.60f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    const AxisSpec axisLambdaMass = {200, 1.f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    const AxisSpec axisVertexZ = {30, -15., 15., "vrtx_{Z} [cm]"};
    const AxisSpec axiscentFT0C = {1020, -1.0, 101.0, "centFT0C(percentile)"};

    const AxisSpec axisP = {200, 0.0f, 10.0f, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPt = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisTPCInnerParam = {200, 0.0f, 10.0f, "#it{p}_{tpcInnerParam} (GeV/#it{c})"};
    const AxisSpec axisTOFExpMom = {200, 0.0f, 10.0f, "#it{p}_{tofExpMom} (GeV/#it{c})"};

    const AxisSpec axisEta = {100, -5, 5, "#eta"};
    const AxisSpec axisPhi = {90, -1, 8, "#phi (radians)"};
    const AxisSpec axisRapidity = {200, -10, 10, "Rapidity (y)"};
    const AxisSpec axisDcaXY = {100, -5, 5, "dcaXY"};
    const AxisSpec axisDcaZ = {100, -5, 5, "dcaZ"};
    const AxisSpec axisDcaXYwide = {2000, -100, 100, "dcaXY"};
    const AxisSpec axisDcaZwide = {2000, -100, 100, "dcaZ"};
    const AxisSpec axisSign = {10, -5, 5, "track.sign"};

    const AxisSpec axisTPCSignal = {100, -1, 1000, "tpcSignal"};
    const AxisSpec axisTOFBeta = {40, -2.0, 2.0, "tofBeta"};

    const AxisSpec axisTPCSignalFine = {10010, -1, 1000, "tpcSignal"};
    const AxisSpec axisTOFBetaFine = {10010, -1, 1000, "tpcSignal"};

    const AxisSpec axisTPCNSigma = {200, -10.0, 10.0, "n#sigma_{TPC}"};
    const AxisSpec axisTOFNSigma = {200, -10.0, 10.0, "n#sigma_{TOF}"};

    const AxisSpec axisTPCNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TPC}^{Pi}"};
    const AxisSpec axisTOFNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TOF}^{Pi}"};

    const AxisSpec axisTPCNClsCrossedRows = {200, -1.5, 198.5, "tpcNClsCrossedRows"};
    const AxisSpec axisIsPVContributor = {4, -1, 3, "isPVContributor"};
    const AxisSpec axisIsGlobalTrack = {4, -1, 3, "isGobalTrack"};
    const AxisSpec axisIsK0sDau = {4, -1, 3, "isK0sDau"};

    const AxisSpec axisDcapostopv = {100000, -50, 50, "dcapostopv"};
    const AxisSpec axisDcanegtopv = {100000, -50, 50, "dcanegtopv"};
    const AxisSpec axisDcaV0daughters = {2000, -10.0, 10.0, "dcaV0daughters"};
    const AxisSpec axisV0cosPA = {3000, -1.5, 1.5, "v0cosPA"};
    const AxisSpec axisV0radius = {100000, -50, 50, "v0radius"};

    const AxisSpec axisParticleCount1 = {60, -10, 50, "particleCount"};
    const AxisSpec axisParticleCount2 = {260, -10, 250, "particleCount"};
    const AxisSpec axisParticleCount3 = {1060, -10, 1050, "particleCount"};

    HistogramConfigSpec histPDcaXY({HistType::kTH2F, {axisP, axisDcaXY}});
    HistogramConfigSpec histPtDcaXY({HistType::kTH2F, {axisPt, axisDcaXY}});
    HistogramConfigSpec histTpcInnerParamDcaXY({HistType::kTH2F, {axisTPCInnerParam, axisDcaXY}});
    HistogramConfigSpec histTofExpMomDcaXY({HistType::kTH2F, {axisTOFExpMom, axisDcaXY}});

    HistogramConfigSpec histPDcaZ({HistType::kTH2F, {axisP, axisDcaZ}});
    HistogramConfigSpec histPtDcaZ({HistType::kTH2F, {axisPt, axisDcaZ}});
    HistogramConfigSpec histTpcInnerParamDcaZ({HistType::kTH2F, {axisTPCInnerParam, axisDcaZ}});
    HistogramConfigSpec histTofExpMomDcaZ({HistType::kTH2F, {axisTOFExpMom, axisDcaZ}});

    HistogramConfigSpec histPPt({HistType::kTH2F, {axisP, axisPt}});
    HistogramConfigSpec histPTpcInnerParam({HistType::kTH2F, {axisP, axisTPCInnerParam}});
    HistogramConfigSpec histPTofExpMom({HistType::kTH2F, {axisP, axisTOFExpMom}});

    HistogramConfigSpec histPTpcSignal({HistType::kTH2F, {axisP, axisTPCSignal}});
    HistogramConfigSpec histTpcInnerParamTpcSignal({HistType::kTH2F, {axisTPCInnerParam, axisTPCSignal}});
    HistogramConfigSpec histTofExpMomTpcSignal({HistType::kTH2F, {axisTOFExpMom, axisTPCSignal}});

    HistogramConfigSpec histPBeta({HistType::kTH2F, {axisP, axisTOFBeta}});
    HistogramConfigSpec histTpcInnerParamBeta({HistType::kTH2F, {axisTPCInnerParam, axisTOFBeta}});
    HistogramConfigSpec histTofExpMomBeta({HistType::kTH2F, {axisTOFExpMom, axisTOFBeta}});

    HistogramConfigSpec histPTpcNSigma({HistType::kTH2F, {axisP, axisTPCNSigma}});
    HistogramConfigSpec histPtTpcNSigma({HistType::kTH2F, {axisPt, axisTPCNSigma}});
    HistogramConfigSpec histTpcInnerParamTpcNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTPCNSigma}});
    HistogramConfigSpec histTofExpMomTpcNSigma({HistType::kTH2F, {axisTOFExpMom, axisTPCNSigma}});
    HistogramConfigSpec histPTofNSigma({HistType::kTH2F, {axisP, axisTOFNSigma}});
    HistogramConfigSpec histPtTofNSigma({HistType::kTH2F, {axisPt, axisTOFNSigma}});
    HistogramConfigSpec histTpcInnerParamTofNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTOFNSigma}});
    HistogramConfigSpec histTofExpMomTofNSigma({HistType::kTH2F, {axisTOFExpMom, axisTOFNSigma}});
    HistogramConfigSpec histTpcNSigmaTofNSigma({HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});

    recoV0s.add("v0Table/Full/h01_K0s_Mass", "K0s_Mass", {HistType::kTH1F, {axisK0sMass}});
    recoV0s.add("v0Table/Full/h02_Lambda_Mass", "Lambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoV0s.add("v0Table/Full/h03_AntiLambda_Mass", "AntiLambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoV0s.add("v0Table/Full/h04_v0DaughterCollisionIndexTag", "hV0s_K0s_v0DaughterCollisionIndexTag", {HistType::kTH1D, {{22, -1.0, 10.0}}});
    recoV0s.add("v0Table/Full/h05_V0Tag", "V0Tag", {HistType::kTH1F, {{12, -2, 10}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda

    // Topological Cuts
    recoV0s.add("v0Table/Full/h06_dcapostopv", "dcapostopv", kTH1F, {axisDcapostopv});
    recoV0s.add("v0Table/Full/h07_dcanegtopv", "dcanegtopv", kTH1F, {axisDcanegtopv});
    recoV0s.add("v0Table/Full/h08_dcaV0daughters", "dcaV0daughters", kTH1F, {axisDcaV0daughters});
    recoV0s.add("v0Table/Full/h09_v0cosPA", "v0cosPA", kTH1F, {axisV0cosPA});
    recoV0s.add("v0Table/Full/h10_v0radius", "v0radius", kTH1F, {axisV0radius});

    // K0s-FullInformation
    recoV0s.add("v0Table/Full/h11_mass", "mass", kTH1F, {axisK0sMass});
    recoV0s.add("v0Table/Full/h12_p", "p", kTH1F, {axisP});
    recoV0s.add("v0Table/Full/h13_pt", "pt", kTH1F, {axisPt});
    recoV0s.add("v0Table/Full/h14_eta", "eta", kTH1F, {axisEta});
    recoV0s.add("v0Table/Full/h15_phi", "phi", kTH1F, {axisPhi});
    recoV0s.add("v0Table/Full/h16_rapidity", "rapidity", kTH1F, {axisRapidity});

    // K0s-Daughter Info
    recoV0s.add("v0Table/Full/Pi/tpcId/h01_p", "p", kTH1F, {axisP});
    recoV0s.add("v0Table/Full/Pi/tpcId/h02_pt", "pt", kTH1F, {axisPt});
    recoV0s.add("v0Table/Full/Pi/tpcId/h03_tpcInnerParam", "tpcInnerParam", kTH1F, {axisTPCInnerParam});
    recoV0s.add("v0Table/Full/Pi/tpcId/h04_tofExpMom", "tofExpMom", kTH1F, {axisTOFExpMom});
    recoV0s.add("v0Table/Full/Pi/tpcId/h05_eta", "eta", kTH1F, {axisEta});
    recoV0s.add("v0Table/Full/Pi/tpcId/h06_phi", "phi", kTH1F, {axisPhi});
    recoV0s.add("v0Table/Full/Pi/tpcId/h07_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoV0s.add("v0Table/Full/Pi/tpcId/h08_isPVContributor", "isPVContributor", kTH1F, {axisIsPVContributor});
    recoV0s.add("v0Table/Full/Pi/tpcId/h09_isGlobalTrack", "isGlobalTrack", kTH1F, {axisIsGlobalTrack});
    recoV0s.add("v0Table/Full/Pi/tpcId/h10_dcaXY", "dcaXY", kTH1F, {axisDcaXY});
    recoV0s.add("v0Table/Full/Pi/tpcId/h11_dcaZ", "dcaZ", kTH1F, {axisDcaZ});

    recoV0s.add("v0Table/Full/Pi/tpcId/h12_p_dcaXY", "p_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoV0s.add("v0Table/Full/Pi/tpcId/h13_p_dcaZ", "p_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoV0s.add("v0Table/Full/Pi/tpcId/h14_pt_dcaXY", "pt_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoV0s.add("v0Table/Full/Pi/tpcId/h15_pt_dcaZ", "pt_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoV0s.add("v0Table/Full/Pi/tpcId/h16_dcaXYwide", "dcaXYwide", kTH1F, {axisDcaXYwide});
    recoV0s.add("v0Table/Full/Pi/tpcId/h17_dcaZwide", "dcaZwide", kTH1F, {axisDcaZwide});

    // K0s-Daughter identification
    // momemtum
    recoV0s.add("v0Table/Full/Pi/tpcId/h20_p_pt", "p_pt", histPPt);
    recoV0s.add("v0Table/Full/Pi/tpcId/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoV0s.add("v0Table/Full/Pi/tpcId/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);
    // tpcSignal
    recoV0s.add("v0Table/Full/Pi/tpcId/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoV0s.add("v0Table/Full/Pi/tpcId/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoV0s.add("v0Table/Full/Pi/tpcId/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);
    // tofBeta
    recoV0s.add("v0Table/Full/Pi/tpcId/h26_p_beta", "p_beta", histPBeta);
    recoV0s.add("v0Table/Full/Pi/tpcId/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoV0s.add("v0Table/Full/Pi/tpcId/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);
    // Look at Pion
    recoV0s.add("v0Table/Full/Pi/tpcId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    recoV0s.addClone("v0Table/Full/Pi/tpcId/", "v0Table/Full/Pi/tofId/"); // for identification using tof+tpc
    recoV0s.addClone("v0Table/Full/Pi/tpcId/", "v0Table/Full/Pi/NoId/");  // for unidentified case // to observe and debug

    recoV0s.addClone("v0Table/Full/", "v0Table/postK0sCheck/");
    recoV0s.addClone("v0Table/Full/", "v0Table/postMassCut/");
    recoV0s.addClone("v0Table/Full/", "v0Table/postSelectionCut/");

    recoV0s.add("v0Table/postSelectionCut/hTrueV0TagCount", "hTrueV0TagCount", {HistType::kTH1F, {{12, -2, 10}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda
    recoV0s.add("v0Table/postSelectionCut/nCommonPionOfDifferentK0s", "nCommonPionOfDifferentK0s", {HistType::kTH1D, {{44, -2, 20}}});

    // Event Selection
    recoEvent.add("recoEvent/h01_CollisionCount", "CollisionCount", {HistType::kTH1D, {{1, 0, 1}}});
    recoEvent.add("recoEvent/h02_VertexXRec", "VertexXRec", {HistType::kTH1D, {{1000, -0.2, 0.2}}});
    recoEvent.add("recoEvent/h03_VertexYRec", "VertexYRec", {HistType::kTH1D, {{1000, -0.2, 0.2}}});
    recoEvent.add("recoEvent/h04_VertexZRec", "VertexZRec", {HistType::kTH1F, {axisVertexZ}});
    recoEvent.add("recoEvent/h05_Centrality", "Centrality", {HistType::kTH1F, {axiscentFT0C}});
    recoEvent.add("recoEvent/h06_V0Size", "V0Size", {HistType::kTH1F, {{60, -10, 50}}});
    recoEvent.add("recoEvent/h07_TracksSize", "TracksSize", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h08_nTrack", "nTrack", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h09_nK0s", "nK0s", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h10_nPiPlus", "nPiPlus", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h11_nPiMinus", "nPiMinus", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h12_nKaPlus", "nKaPlus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h13_nKaMinus", "nKaMinus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h14_nProton", "nProton", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h15_nPBar", "nPBar", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h16_nElPlus", "nElPlus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h17_nElMinus", "nElMinus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h18_nDePlus", "nDePlus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h19_nDeMinus", "nDeMinus", {HistType::kTH1F, {axisParticleCount1}});

    //
    // K0s reconstruction
    recoK0s.add("recoK0s/PreSel/h01_K0s_Mass", "K0s_Mass", {HistType::kTH1F, {axisK0sMass}});
    recoK0s.add("recoK0s/PreSel/h02_Lambda_Mass", "Lambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoK0s.add("recoK0s/PreSel/h03_AntiLambda_Mass", "AntiLambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoK0s.add("recoK0s/PreSel/h04_v0DaughterCollisionIndexTag", "hV0s_K0s_v0DaughterCollisionIndexTag", {HistType::kTH1D, {{22, -1.0, 10.0}}});
    recoK0s.add("recoK0s/PreSel/h05_V0Tag", "V0Tag", {HistType::kTH1F, {{12, -2, 10}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda

    // Topological Cuts
    recoK0s.add("recoK0s/PreSel/h06_dcapostopv", "dcapostopv", kTH1F, {axisDcapostopv});
    recoK0s.add("recoK0s/PreSel/h07_dcanegtopv", "dcanegtopv", kTH1F, {axisDcanegtopv});
    recoK0s.add("recoK0s/PreSel/h08_dcaV0daughters", "dcaV0daughters", kTH1F, {axisDcaV0daughters});
    recoK0s.add("recoK0s/PreSel/h09_v0cosPA", "v0cosPA", kTH1F, {axisV0cosPA});
    recoK0s.add("recoK0s/PreSel/h10_v0radius", "v0radius", kTH1F, {axisV0radius});

    // K0s-FullInformation
    recoK0s.add("recoK0s/PreSel/h11_mass", "mass", kTH1F, {axisK0sMass});
    recoK0s.add("recoK0s/PreSel/h12_p", "p", kTH1F, {axisP});
    recoK0s.add("recoK0s/PreSel/h13_pt", "pt", kTH1F, {axisPt});
    recoK0s.add("recoK0s/PreSel/h14_eta", "eta", kTH1F, {axisEta});
    recoK0s.add("recoK0s/PreSel/h15_phi", "phi", kTH1F, {axisPhi});
    recoK0s.add("recoK0s/PreSel/h16_rapidity", "rapidity", kTH1F, {axisRapidity});

    // K0s-Daughter Info
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h01_p", "p", kTH1F, {axisP});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h02_pt", "pt", kTH1F, {axisPt});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h03_tpcInnerParam", "tpcInnerParam", kTH1F, {axisTPCInnerParam});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h04_tofExpMom", "tofExpMom", kTH1F, {axisTOFExpMom});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h05_eta", "eta", kTH1F, {axisEta});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h06_phi", "phi", kTH1F, {axisPhi});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h07_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h08_isPVContributor", "isPVContributor", kTH1F, {axisIsPVContributor});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h09_isGlobalTrack", "isGlobalTrack", kTH1F, {axisIsGlobalTrack});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h10_dcaXY", "dcaXY", kTH1F, {axisDcaXY});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h11_dcaZ", "dcaZ", kTH1F, {axisDcaZ});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h12_p_dcaXY", "p_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h13_p_dcaZ", "p_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h14_pt_dcaXY", "pt_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h15_pt_dcaZ", "pt_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h16_dcaXYwide", "dcaXYwide", kTH1F, {axisDcaXYwide});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h17_dcaZwide", "dcaZwide", kTH1F, {axisDcaZwide});

    // K0s-Daughter identification
    // momemtum
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h20_p_pt", "p_pt", histPPt);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);
    // tpcSignal
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);
    // tofBeta
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h26_p_beta", "p_beta", histPBeta);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);
    // Look at Pion
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    recoK0s.addClone("recoK0s/PreSel/Pi/tpcId/", "recoK0s/PreSel/Pi/tofId/"); // for identification using tof+tpc
    recoK0s.addClone("recoK0s/PreSel/Pi/tpcId/", "recoK0s/PreSel/Pi/NoId/");  // for unidentified case // to observe and debug

    recoK0s.addClone("recoK0s/PreSel/", "recoK0s/PostSel/"); // for unidentified case // to observe and debug

    recoK0s.add("recoK0s/PostSel/mK0s_vs_centFTOC", "mK0s_vs_centFTOC", kTH2F, {axiscentFT0C, axisK0sMass});

    // Tracks reconstruction
    // FullTrack
    recoTracks.add("recoTracks/PreSel/h01_p", "p", {HistType::kTH1F, {axisP}});
    recoTracks.add("recoTracks/PreSel/h02_pt", "pt", {HistType::kTH1F, {axisPt}});
    recoTracks.add("recoTracks/PreSel/h03_tpcInnerParam", "tpcInnerParam", {HistType::kTH1F, {axisTPCInnerParam}});
    recoTracks.add("recoTracks/PreSel/h04_tofExpMom", "tofExpMom", {HistType::kTH1F, {axisTOFExpMom}});
    recoTracks.add("recoTracks/PreSel/h05_eta", "eta", {HistType::kTH1F, {axisEta}});
    recoTracks.add("recoTracks/PreSel/h06_phi", "phi", {HistType::kTH1F, {axisPhi}});
    recoTracks.add("recoTracks/PreSel/h07_dcaXY", "dcaXY", {HistType::kTH1F, {axisDcaXY}});
    recoTracks.add("recoTracks/PreSel/h08_dcaZ", "dcaZ", {HistType::kTH1F, {axisDcaZ}});
    recoTracks.add("recoTracks/PreSel/h09_sign", "sign", {HistType::kTH1D, {axisSign}});

    // DcaXY
    recoTracks.add("recoTracks/PreSel/h10_p_dcaXY", "p_dcaXY", histPDcaXY);
    recoTracks.add("recoTracks/PreSel/h11_pt_dcaXY", "pt_dcaXY", histPtDcaXY);
    recoTracks.add("recoTracks/PreSel/h12_tpcInnerParam_dcaXY", "tpcInnerParam_dcaXY", histTpcInnerParamDcaXY);
    recoTracks.add("recoTracks/PreSel/h13_tofExpMom_dcaXY", "tofExpMom_dcaXY", histTofExpMomDcaXY);

    // DcaZ
    recoTracks.add("recoTracks/PreSel/h14_p_dcaZ", "p_dcaZ", histPDcaZ);
    recoTracks.add("recoTracks/PreSel/h15_pt_dcaZ", "pt_dcaZ", histPtDcaZ);
    recoTracks.add("recoTracks/PreSel/h16_tpcInnerParam_dcaZ", "tpcInnerParam_dcaZ", histTpcInnerParamDcaZ);
    recoTracks.add("recoTracks/PreSel/h17_tofExpMom_dcaZ", "tofExpMom_dcaZ", histTofExpMomDcaZ);

    // momemtum
    recoTracks.add("recoTracks/PreSel/h20_p_pt", "p_pt", histPPt);
    recoTracks.add("recoTracks/PreSel/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoTracks.add("recoTracks/PreSel/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);

    // tpcSignal
    recoTracks.add("recoTracks/PreSel/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoTracks.add("recoTracks/PreSel/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoTracks.add("recoTracks/PreSel/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);

    // tofBeta
    recoTracks.add("recoTracks/PreSel/h26_p_beta", "p_beta", histPBeta);
    recoTracks.add("recoTracks/PreSel/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoTracks.add("recoTracks/PreSel/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);

    // Look at Pion
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);
    // Pion

    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/Ka/"); // Kaon
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/Pr/"); // Proton
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/El/"); // Electron
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/De/"); // Deuteron

    // Write Code for naming the axis for Identified Particles

    // Selection
    recoTracks.addClone("recoTracks/PreSel/", "recoTracks/PostSel/");
    //

    // Analysis
    // momemtum
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h20_p_pt", "p_pt", histPPt);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);
    // tpcSignal
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);
    // tofBeta
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h26_p_beta", "p_beta", histPBeta);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);
    // Pion
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    recoAnalysis.addClone("recoAnalysis/Pi/tpcId/", "recoAnalysis/Pi/tofId/");
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Ka/"); // Kaon
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Pr/"); // Proton
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/El/"); // Electron
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/De/"); // Deuteron

    recoAnalysis.add("recoAnalysis/SelectedTrack_IdentificationTag", "SelectedTrack_IdentificationTag", kTH1D, {{34, -1.5, 32.5, "trackTAG"}});
    recoAnalysis.add("recoAnalysis/RejectedTrack_RejectionTag", "RejectedTrack_RejectionTag", kTH1D, {{16, -1.5, 6.5, "rejectionTAG"}});

    recoAnalysis.add("recoAnalysis/Sparse_Full_K0sPiKa", "Sparse_Full_K0sPiKa", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_K0sPrDe", "Sparse_Full_K0sPrDe", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_K0sKaEl", "Sparse_Full_K0sKaEl", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_PiKaPr", "Sparse_Full_PiKaPr", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_PiElDe", "Sparse_Full_PiElDe", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_KaPrDe", "Sparse_Full_KaPrDe", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_PrElDe", "Sparse_Full_PrElDe", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});

    recoAnalysis.add("recoAnalysis/Sparse_newDynm_K0s_Ka", "Sparse_newDynm_K0s_Ka", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {500, -1.5, 498.5, "nKaon"}, {10000, -1.5, 9998.5, "(nK0s)^{2}"}, {250000, -1.5, 249998.5, "(nKaon)^{2}"}, {500, -1.5, 498.5, "(nK0s*nKaon)"}});
    recoAnalysis.add("recoAnalysis/Sparse_newDynm_Kp_Km", "Sparse_newDynm_Kp_Km", kTHnSparseD, {axiscentFT0C, {2000, -1.5, 1998.5, "nTrack"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {250000, -1.5, 249998.5, "(nKaPlus)^{2}"}, {250000, -1.5, 249998.5, "(nKaMinus)^{2}"}, {250000, -1.5, 249998.5, "(nKaPlus*nKaMinus)"}});
    //

    // Printing the Stored Registry information
    LOG(info) << "Printing Stored Registry Information";
    LOG(info) << " DEBUG :: 01- recoV0s.print()";
    recoV0s.print();
    LOG(info) << " DEBUG :: 02- recoEvent.print()";
    recoEvent.print();
    LOG(info) << " DEBUG :: 03- recoK0s.print()";
    recoK0s.print();
    LOG(info) << " DEBUG :: 04- recoTracks.print()";
    recoTracks.print();
    LOG(info) << " DEBUG :: 05- recoAnalysis.print()";
    recoAnalysis.print();
  }

  // tpc Selections
  template <typename T>
  bool selPionTPCInnerParam(T track)
  {
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaKa()) > 3.0 && std::abs(track.tpcNSigmaPr()) > 3.0 && std::abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && std::abs(track.tpcNSigmaPi()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPi()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selKaonTPCInnerParam(T track)
  {
    // p dependent cuts
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaPi()) > 3.0 && std::abs(track.tpcNSigmaPr()) > 3.0 && std::abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && std::abs(track.tpcNSigmaKa()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaKa()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selProtonTPCInnerParam(T track)
  {
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaPi()) > 3.0 && std::abs(track.tpcNSigmaKa()) > 3.0 && std::abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.60 && std::abs(track.tpcNSigmaPr()) < 3.0) {
        return true;
      }
      if (1.60 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPr()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selDeuteronTPCInnerParam(T track)
  {
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaPi()) > 3.0 && std::abs(track.tpcNSigmaKa()) > 3.0 && std::abs(track.tpcNSigmaPr()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.80 && std::abs(track.tpcNSigmaDe()) < 3.0) {
        return true;
      }
      if (1.80 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaDe()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selElectronTPCInnerParam(T track)
  {
    if (track.tpcNSigmaEl() < 3.0 && track.tpcNSigmaPi() > 3.0 && track.tpcNSigmaKa() > 3.0 && track.tpcNSigmaPr() > 3.0 && track.tpcNSigmaDe() > 3.0) {
      return true;
    }
    return false;
  }
  //

  // TOF Selections
  // Pion
  template <typename T>
  bool selPionTOF(T track)
  {
    if (track.p() <= 0.75 && std::abs(track.tpcNSigmaPi()) < 3.0 && std::abs(track.tofNSigmaPi()) < 3.0 && std::abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    } else if (0.75 < track.p() // after p = 0.75, Pi and Ka lines of nSigma 3.0 will start intersecting
               && std::abs(track.tpcNSigmaPi()) < 2.0 && std::abs(track.tofNSigmaPi()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaonTOF(T track)
  {
    if (track.p() <= 0.75 && std::abs(track.tpcNSigmaKa()) < 3.0 && std::abs(track.tofNSigmaKa()) < 3.0) {
      return true;
    }
    if (0.75 < track.p() && track.p() <= 1.30 // after 0.75 Pi and Ka lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaKa()) < 3.0 && std::abs(track.tofNSigmaKa()) < 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    }
    if (1.30 < track.p() // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaKa()) < 2.0 && std::abs(track.tofNSigmaKa()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProtonTOF(T track)
  {
    if (track.p() <= 1.30 && std::abs(track.tpcNSigmaPr()) < 3.0 && std::abs(track.tofNSigmaPr()) < 3.0) {
      return true;
    }
    if (1.30 < track.p() && track.p() <= 3.10                                                                                                                                                                                                     // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaPr()) < 3.0 && std::abs(track.tofNSigmaPr()) < 3.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0 // Some Deuteron contamination is still coming in p dependent cuts
    ) {
      return true;
    }
    if (3.10 < track.p() // after 3.10 Pr and De lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaPr()) < 2.0 && std::abs(track.tofNSigmaPr()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteronTOF(T track)
  {
    if (track.p() <= 3.10 && std::abs(track.tpcNSigmaDe()) < 3.0 && std::abs(track.tofNSigmaDe()) < 3.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0) {
      return true;
    }
    if (3.10 < track.p() // after 3.10 De and Pr lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaDe()) < 2.0 && std::abs(track.tofNSigmaDe()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0) {
      return true;
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectronTOF(T track)
  {
    if (
      (std::pow(track.tpcNSigmaEl(), 2) + std::pow(track.tofNSigmaEl(), 2)) < 9.00 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }
  //

  // SelectionFunctions
  // Pion
  template <typename T>
  bool selPion(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selPionTPCInnerParam(track)) {
      IdMethod = 1;
      return selPionTOF(track);
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaon(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selKaonTPCInnerParam(track)) {
      IdMethod = 1;
      return selKaonTOF(track);
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProton(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selProtonTPCInnerParam(track)) {
      IdMethod = 1;
      return selProtonTOF(track);
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteron(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selDeuteronTPCInnerParam(track)) {
      IdMethod = 1;
      return selDeuteronTOF(track);
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectron(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selElectronTPCInnerParam(track)) {
      IdMethod = 1;
      return selElectronTOF(track);
    }
    return false;
  }

  template <typename T>
  int findV0Tag(const T& posDaughterTrack, const T& negDaughterTrack, int& posPiIdMethod, int& posPrIdMethod, int& negPiIdMethod, int& negPrIdMethod)
  {
    bool posIsPion = false;
    bool posIsProton = false;
    bool negIsPion = false;
    bool negIsProton = false;

    int v0TagValue = 0;

    // Check if positive track is pion or proton
    if (selPion(posDaughterTrack, posPiIdMethod))
      posIsPion = true; // Coming From K0s    -> PiPlus + PiMinus and AntiLambda -> PiPlus + AntiProton
    if (selProton(posDaughterTrack, posPrIdMethod))
      posIsProton = true; // Coming From Lambda -> proton + PiMinus
    if (selPion(negDaughterTrack, negPiIdMethod))
      negIsPion = true; // Coming From K0s       -> PiPlus + PiMinus and Lambda -> proton + PiMinus
    if (selProton(negDaughterTrack, negPrIdMethod))
      negIsProton = true; // Coming From AntiLambda -> PiPlus + AntiProton

    if (posIsPion && negIsPion) {
      v0TagValue += 1; // It is K0s
    }
    if (posIsProton && negIsPion) {
      v0TagValue += 2; // It is Lambda
    }
    if (posIsPion && negIsProton) {
      v0TagValue += 4; // It is AntiLambda
    }
    return v0TagValue;
  }

  template <typename T, typename U>
  int findCollisionIndexTag(const T& v0, const U& posDaughterTrack, const U& negDaughterTrack)
  {
    int v0daughterCollisionIndexTag = 0;
    if (v0.collisionId() == posDaughterTrack.collisionId()) {
      v0daughterCollisionIndexTag = +1;
    }
    if (v0.collisionId() == negDaughterTrack.collisionId()) {
      v0daughterCollisionIndexTag = +2;
    }
    if (posDaughterTrack.collisionId() == negDaughterTrack.collisionId()) {
      v0daughterCollisionIndexTag = +4;
    }
    return v0daughterCollisionIndexTag;
  }

  template <typename T>
  bool selK0s(T v0)
  {
    if (mLowK0s < v0.mK0Short() && v0.mK0Short() < mHighK0s && 0.1 < v0.pt() && v0.pt() < 1.5 && std::abs(v0.rapidity(MassK0Short)) < 0.5) {
      return true;
    } else {
      return false;
    }
  }

  template <typename T>
  void findRepeatedEntries(std::vector<int64_t> ParticleList, T hist)
  {
    for (uint ii = 0; ii < ParticleList.size(); ii++) {
      int nCommonCount = 0; // checking the repeat number of track
      for (uint jj = 0; jj < ParticleList.size(); jj++) {
        if (ParticleList[jj] == ParticleList[ii]) {
          if (jj < ii) {
            break;
          } // break if it was already counted
          nCommonCount++; // To Calculate no of times the entry was repeated
        }
      }
      hist->Fill(nCommonCount);
    }
  }

  template <typename T>
  bool checkTrackSelection(const T& track, const std::vector<int64_t>& DauParticleList, uint& skippingPosition, int& rejectionTag)
  {
    if (track.tpcNClsCrossedRows() < 70) {
      rejectionTag = 1;
      return false;
    }
    if (std::fabs(track.dcaXY()) > 0.2) {
      rejectionTag = 2;
      return false;
    }
    if (!track.isGlobalTrack()) {
      rejectionTag = 3;
      return false;
    }

    bool flagDaughterTrack = false;
    if (track.globalIndex() == DauParticleList[skippingPosition]) {
      flagDaughterTrack = true;
      skippingPosition++;
    }
    if (flagDaughterTrack) {
      rejectionTag = 4;
      return false;
    }
    // Event and Pt filter will filter out some tracks and collisions used for v0 reconstruction.
    if (track.globalIndex() > DauParticleList[skippingPosition] && skippingPosition < DauParticleList.size()) {
      skippingPosition++;
    }
    return true;
  }

  enum HistRegEnum {
    v0TableFull = 0,
    v0TablePostK0sCheck,
    v0TablePostMassCut,
    v0TablePostSelectionCut,
    recoK0sPreSel,
    recoK0sPostSel,
    recoTrackPreSel,
    recoTrackPostSel,
    recoAnalysisDir
  };

  static constexpr std::string_view HistRegDire[] = {
    "v0Table/Full/",
    "v0Table/postK0sCheck/",
    "v0Table/postMassCut/",
    "v0Table/postSelectionCut/",
    "recoK0s/PreSel/",
    "recoK0s/PostSel/",
    "recoTracks/PreSel/",
    "recoTracks/PostSel/",
    "recoAnalysis/"};

  enum PidEnum {
    kPi = 0, // dont use kPion, kKaon, as these enumeration
    kKa,     // are already defined in $ROOTSYS/root/include/TPDGCode.h
    kPr,
    kEl,
    kDe
  };

  static constexpr std::string_view PidDire[] = {
    "Pi/",
    "Ka/",
    "Pr/",
    "El/",
    "De/"};

  enum DetEnum {
    tpcId = 0,
    tofId,
    NoId
  };

  static constexpr std::string_view DetDire[] = {
    "tpcId/",
    "tofId/",
    "NoId/"};

  template <int Mode, int pidMode, int detMode, bool fillSignal, typename H, typename T>
  void fillIdentificationQA(H histReg, const T& track)
  {
    float tpcNSigmaVal = -999, tofNSigmaVal = -999;
    switch (pidMode) {
      case kPi:
        tpcNSigmaVal = track.tpcNSigmaPi();
        tofNSigmaVal = track.tofNSigmaPi();
        break;
      case kKa:
        tpcNSigmaVal = track.tpcNSigmaKa();
        tofNSigmaVal = track.tofNSigmaKa();
        break;
      case kPr:
        tpcNSigmaVal = track.tpcNSigmaPr();
        tofNSigmaVal = track.tofNSigmaPr();
        break;
      case kEl:
        tpcNSigmaVal = track.tpcNSigmaEl();
        tofNSigmaVal = track.tofNSigmaEl();
        break;
      case kDe:
        tpcNSigmaVal = track.tpcNSigmaDe();
        tofNSigmaVal = track.tofNSigmaDe();
        break;
      default:
        tpcNSigmaVal = -999, tofNSigmaVal = -999;
        break;
    }

    if (fillSignal) {
      // momemtum
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h20_p_pt"), track.p(), track.pt());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h21_p_tpcInnerParam"), track.p(), track.tpcInnerParam());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h22_p_tofExpMom"), track.p(), track.tofExpMom());
      // tpcSignal
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h23_p_tpcSignal"), track.p(), track.tpcSignal());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h24_tpcInnerParam_tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h25_tofExpMom_tpcSignal"), track.tofExpMom(), track.tpcSignal());
      // tofBeta
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h26_p_beta"), track.p(), track.beta());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h27_tpcInnerParam_beta"), track.tpcInnerParam(), track.beta());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h28_tofExpMom_beta"), track.tofExpMom(), track.beta());
    }
    // NSigma
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h29_p_tpcNSigma"), track.p(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h30_pt_tpcNSigma"), track.pt(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h31_tpcInnerParam_tpcNSigma"), track.tpcInnerParam(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h32_tofExpMom_tpcNSigma"), track.tofExpMom(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h33_p_tofNSigma"), track.p(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h34_pt_tofNSigma"), track.pt(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h35_tpcInnerParam_tofNSigma"), track.tpcInnerParam(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h36_tofExpMom_tofNSigma"), track.tofExpMom(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h37_tpcNSigma_tofNSigma"), tpcNSigmaVal, tofNSigmaVal);
  }

  template <int Mode, typename T>
  void fillTrackQA(T track)
  {
    // FullTrack
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h01_p"), track.p());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h02_pt"), track.pt());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h03_tpcInnerParam"), track.tpcInnerParam());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h04_tofExpMom"), track.tofExpMom());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h05_eta"), track.eta());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h06_phi"), track.phi());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h07_dcaXY"), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h08_dcaZ"), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h09_sign"), track.sign());
    // DcaXY
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h10_p_dcaXY"), track.p(), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h11_pt_dcaXY"), track.pt(), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h12_tpcInnerParam_dcaXY"), track.tpcInnerParam(), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h13_tofExpMom_dcaXY"), track.tofExpMom(), track.dcaXY());

    // DcaZ
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h14_p_dcaZ"), track.p(), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h15_pt_dcaZ"), track.pt(), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h16_tpcInnerParam_dcaZ"), track.tpcInnerParam(), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h17_tofExpMom_dcaZ"), track.tofExpMom(), track.dcaZ());

    // momemtum
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h20_p_pt"), track.p(), track.pt());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h21_p_tpcInnerParam"), track.p(), track.tpcInnerParam());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h22_p_tofExpMom"), track.p(), track.tofExpMom());

    // tpcSignal
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h23_p_tpcSignal"), track.p(), track.tpcSignal());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h24_tpcInnerParam_tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h25_tofExpMom_tpcSignal"), track.tofExpMom(), track.tpcSignal());

    // tofBeta
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h26_p_beta"), track.p(), track.beta());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h27_tpcInnerParam_beta"), track.tpcInnerParam(), track.beta());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h28_tofExpMom_beta"), track.tofExpMom(), track.beta());

    fillIdentificationQA<Mode, kPi, NoId, false>(recoTracks, track); // Look at Pion
    fillIdentificationQA<Mode, kKa, NoId, false>(recoTracks, track); // Look at Kaon
    fillIdentificationQA<Mode, kPr, NoId, false>(recoTracks, track); // Look at Proton
    fillIdentificationQA<Mode, kEl, NoId, false>(recoTracks, track); // Look at Electron
    fillIdentificationQA<Mode, kDe, NoId, false>(recoTracks, track); // Look at Deuteron
  }

  template <int Mode, int pidMode, int detMode, typename H, typename T>
  void fillV0DaughterQA(H histReg, const T& track, double particleMass)
  {
    // K0s-Daughter Info
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h01_p"), track.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h02_pt"), track.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h03_tpcInnerParam"), track.tpcInnerParam());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h04_tofExpMom"), track.tofExpMom());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h05_eta"), track.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h06_phi"), track.phi());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h07_rapidity"), track.rapidity(particleMass));
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h08_isPVContributor"), track.isPVContributor());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h09_isGlobalTrack"), track.isGlobalTrack());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h10_dcaXY"), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h11_dcaZ"), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h12_p_dcaXY"), track.p(), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h13_p_dcaZ"), track.p(), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h14_pt_dcaXY"), track.pt(), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h15_pt_dcaZ"), track.pt(), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h16_dcaXYwide"), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h17_dcaZwide"), track.dcaZ());

    fillIdentificationQA<Mode, kPi, detMode, true>(histReg, track);
  }

  template <int Mode, typename H, typename T, typename U>
  void fillV0QA(H histReg, const T& v0, const U& posDaughterTrack, const U& negDaughterTrack, const int& v0Tag, const int& v0DauCollisionIndexTag, const int& posPiIdMethod, const int& negPiIdMethod)
  {
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h01_K0s_Mass"), v0.mK0Short());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h02_Lambda_Mass"), v0.mLambda());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h03_AntiLambda_Mass"), v0.mAntiLambda());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h04_v0DaughterCollisionIndexTag"), v0DauCollisionIndexTag);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h05_V0Tag"), v0Tag);

    // Topological Cuts
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h06_dcapostopv"), v0.dcapostopv());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h07_dcanegtopv"), v0.dcanegtopv());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h08_dcaV0daughters"), v0.dcaV0daughters());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h09_v0cosPA"), v0.v0cosPA());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h10_v0radius"), v0.v0radius());

    // K0s-FullInformation
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mK0Short());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h12_p"), v0.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h13_pt"), v0.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h14_eta"), v0.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h15_phi"), v0.phi());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity"), v0.rapidity(MassK0Short));

    if (posPiIdMethod == 0) {
      fillV0DaughterQA<Mode, kPi, tpcId>(histReg, posDaughterTrack, MassProton);
    } else if (posPiIdMethod == 1) {
      fillV0DaughterQA<Mode, kPi, tofId>(histReg, posDaughterTrack, MassProton);
    } else if (posPiIdMethod == -1) {
      fillV0DaughterQA<Mode, kPi, NoId>(histReg, posDaughterTrack, MassProton);
    }

    if (negPiIdMethod == 0) {
      fillV0DaughterQA<Mode, kPi, tpcId>(histReg, negDaughterTrack, MassProton);
    } else if (negPiIdMethod == 1) {
      fillV0DaughterQA<Mode, kPi, tofId>(histReg, negDaughterTrack, MassProton);
    } else if (negPiIdMethod == -1) {
      fillV0DaughterQA<Mode, kPi, NoId>(histReg, negDaughterTrack, MassProton);
    }
  }

  // Event Filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZvertex);

  // Track Filter
  Filter ptFilter = (o2::aod::track::pt) > 0.15f && (o2::aod::track::pt) < 2.0f;

  // Filters on V0s
  Filter preFilterv0 = (nabs(aod::v0data::dcapostopv) > v0settingDcaPosToPV &&
                        nabs(aod::v0data::dcanegtopv) > v0settingDcaNegToPV &&
                        aod::v0data::dcaV0daughters < v0settingDcaV0Dau);

  using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                               aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;

  using MyTracks = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>>;

  using MyV0s = soa::Filtered<aod::V0Datas>;

  // For manual sliceBy
  Preslice<MyTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<MyV0s> v0sPerCollisionPreslice = o2::aod::track::collisionId;

  // Declaring vectors outside the process to avoid slight overhead for stack allocation and deallocation during each iteration.
  std::vector<int64_t> k0sPosDauList;
  std::vector<int64_t> k0sNegDauList;

  // fullDauList
  std::vector<int64_t> fullDauList;

  void process(MyCollisions const& collisions,
               MyV0s const& V0s,
               MyTracks const& tracks)
  {

    k0sPosDauList.clear();
    k0sNegDauList.clear();
    fullDauList.clear();

    int posPiIdMethod = -1;
    int posPrIdMethod = -1;
    int negPiIdMethod = -1;
    int negPrIdMethod = -1;

    bool isK0s = false;

    int v0Tag = 0;
    int trueV0TagValue = 0;
    int v0DauCollisionIndexTag = 0;

    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<MyTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<MyTracks>();

      posPiIdMethod = -1;
      posPrIdMethod = -1;
      negPiIdMethod = -1;
      negPrIdMethod = -1;
      v0Tag = findV0Tag(posDaughterTrack, negDaughterTrack, posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod);
      v0DauCollisionIndexTag = findCollisionIndexTag(v0, posDaughterTrack, negDaughterTrack);

      isK0s = false;
      if ((v0Tag & 1) == 1)
        isK0s = true;
      trueV0TagValue = 0;

      fillV0QA<v0TableFull>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);

      // cut on dynamic columns for v0 particles
      if (v0.v0cosPA() < v0settingCosPA)
        continue;
      if (v0.v0radius() < v0settingRadius)
        continue;

      // K0s Analysis
      if (isK0s) {
        fillV0QA<v0TablePostK0sCheck>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
        // K0s mass cut
        if (mLowK0s < v0.mK0Short() && v0.mK0Short() < mHighK0s) {
          fillV0QA<v0TablePostMassCut>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
        }

        // Final K0s Selection.
        if (selK0s(v0)) {
          fillV0QA<v0TablePostSelectionCut>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
          trueV0TagValue += 1;
          k0sPosDauList.push_back(posDaughterTrack.globalIndex());
          k0sNegDauList.push_back(negDaughterTrack.globalIndex());
        }
        recoV0s.fill(HIST(HistRegDire[v0TablePostSelectionCut]) + HIST("hTrueV0TagCount"), trueV0TagValue); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda
      } // End of K0s block
    } // End of V0s Loop

    findRepeatedEntries(k0sPosDauList, recoV0s.get<TH1>(HIST(HistRegDire[v0TablePostSelectionCut]) + HIST("nCommonPionOfDifferentK0s")));
    findRepeatedEntries(k0sNegDauList, recoV0s.get<TH1>(HIST(HistRegDire[v0TablePostSelectionCut]) + HIST("nCommonPionOfDifferentK0s")));

    // Obtain one single new daughter vector to remove double counting
    fullDauList.insert(fullDauList.end(), k0sPosDauList.begin(), k0sPosDauList.end());
    fullDauList.insert(fullDauList.end(), k0sNegDauList.begin(), k0sNegDauList.end());

    // Sort and Remove repeated entries
    std::sort(fullDauList.begin(), fullDauList.end());
    auto last = std::unique(fullDauList.begin(), fullDauList.end()); // std::unique only moves duplicates to end of the vector
    fullDauList.erase(last, fullDauList.end());                      // last is the iterator position from where duplicate entries start

    // Check sorting
    if (!std::is_sorted(fullDauList.begin(), fullDauList.end())) {
      LOG(error) << "fullDauList is unsorted, will give wrong results when v0 and collisions will be checked";
    }

    // Declaring variables outside the loop to avoid slight overhead for stack allocation and deallocation during each iteration.
    uint skippingPosition = 0;
    int nK0s = 0;
    int nPiPlus = 0;
    int nPiMinus = 0;
    int nKaPlus = 0;
    int nKaMinus = 0;
    int nProton = 0;
    int nPBar = 0;
    int nElPlus = 0;
    int nElMinus = 0;
    int nDePlus = 0;
    int nDeMinus = 0;

    int nTrack = 0;
    int nKaon = 0;
    double centFT0C = 0;

    int nRejectedPiPlus = 0;
    int nRejectedPiMinus = 0;
    int rejectionTag = 0;

    bool trackIsPion = false;
    bool trackIsKaon = false;
    bool trackIsProton = false;
    bool trackIsElectron = false;
    bool trackIsDeuteron = false;

    int trackIdTag = 0;
    int idMethodPi = -1;
    int idMethodKa = -1;
    int idMethodPr = -1;
    int idMethodEl = -1;
    int idMethodDe = -1;

    for (const auto& collision : collisions) {

      nK0s = 0;
      nPiPlus = 0;
      nPiMinus = 0;
      nKaPlus = 0;
      nKaMinus = 0;
      nProton = 0;
      nPBar = 0;
      nElPlus = 0;
      nElMinus = 0;
      nDePlus = 0;
      nDeMinus = 0;
      nTrack = 0;
      nKaon = 0;

      centFT0C = collision.centFT0C();

      // group tracks, v0s manually
      const uint64_t collIdx = collision.globalIndex();
      const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);
      const auto v0sTablePerColl = V0s.sliceBy(v0sPerCollisionPreslice, collIdx);

      for (const auto& v0 : v0sTablePerColl) {
        if (v0.v0cosPA() < v0settingCosPA)
          continue; // cut on dynamic columns for v0 particles
        if (v0.v0radius() < v0settingRadius)
          continue;

        const auto& posDaughterTrack = v0.posTrack_as<MyTracks>();
        const auto& negDaughterTrack = v0.negTrack_as<MyTracks>();

        bool isK0s = false;

        posPiIdMethod = -1;
        posPrIdMethod = -1;
        negPiIdMethod = -1;
        negPrIdMethod = -1;
        v0Tag = findV0Tag(posDaughterTrack, negDaughterTrack, posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod);
        v0DauCollisionIndexTag = findCollisionIndexTag(v0, posDaughterTrack, negDaughterTrack);

        if ((v0Tag & 1) == 1)
          isK0s = true;

        fillV0QA<recoK0sPreSel>(recoK0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);

        // K0s Analysis
        if (isK0s && selK0s(v0)) {
          fillV0QA<recoK0sPostSel>(recoK0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
          recoK0s.fill(HIST(HistRegDire[recoK0sPostSel]) + HIST("mK0s_vs_centFTOC"), collision.centFT0C(), v0.mK0Short()); // centrality dependent mass
          nK0s++;
        } // End of K0s block
      } // End of V0s Loop

      nTrack = 0;
      nRejectedPiPlus = 0;
      nRejectedPiMinus = 0;
      for (const auto& track : tracksTablePerColl) {

        fillTrackQA<recoTrackPreSel>(track);
        rejectionTag = 0;
        if (!checkTrackSelection(track, fullDauList, skippingPosition, rejectionTag)) {
          if (rejectionTag == 4) {
            if (track.sign() > 0) {
              nRejectedPiPlus++;
            }
            if (track.sign() < 0) {
              nRejectedPiMinus++;
            }
          }
          recoAnalysis.fill(HIST("recoAnalysis/RejectedTrack_RejectionTag"), rejectionTag);
          continue;
        }

        fillTrackQA<recoTrackPostSel>(track);

        nTrack++;
        // Do Proper Track Identification
        trackIsPion = false;
        trackIsKaon = false;
        trackIsProton = false;
        trackIsElectron = false;
        trackIsDeuteron = false;

        trackIdTag = 0;
        idMethodPi = -1;
        idMethodKa = -1;
        idMethodPr = -1;
        idMethodEl = -1;
        idMethodDe = -1;

        if (selPion(track, idMethodPi)) {
          trackIsPion = true;
          trackIdTag += 1;
        }
        if (selKaon(track, idMethodKa)) {
          trackIsKaon = true;
          trackIdTag += 2;
        }
        if (selProton(track, idMethodPr)) {
          trackIsProton = true;
          trackIdTag += 4;
        }
        if (selElectron(track, idMethodEl)) {
          trackIsElectron = true;
          trackIdTag += 8;
        }
        if (selDeuteron(track, idMethodDe)) {
          trackIsDeuteron = true;
          trackIdTag += 16;
        }

        if (trackIsPion) {
          if (idMethodPi == 0)
            fillIdentificationQA<recoAnalysisDir, kPi, tpcId, true>(recoAnalysis, track);
          if (idMethodPi == 1)
            fillIdentificationQA<recoAnalysisDir, kPi, tofId, true>(recoAnalysis, track);
          if (track.sign() > 0) {
            nPiPlus++;
          }
          if (track.sign() < 0) {
            nPiMinus++;
          }
        }
        if (trackIsKaon) {
          if (idMethodKa == 0)
            fillIdentificationQA<recoAnalysisDir, kKa, tpcId, true>(recoAnalysis, track);
          if (idMethodKa == 1)
            fillIdentificationQA<recoAnalysisDir, kKa, tofId, true>(recoAnalysis, track);

          if (track.sign() > 0) {
            nKaPlus++;
          }
          if (track.sign() < 0) {
            nKaMinus++;
          }
        }
        if (trackIsProton) {
          if (idMethodPr == 0)
            fillIdentificationQA<recoAnalysisDir, kPr, tpcId, true>(recoAnalysis, track);
          if (idMethodPr == 1)
            fillIdentificationQA<recoAnalysisDir, kPr, tofId, true>(recoAnalysis, track);

          if (track.sign() > 0) {
            nProton++;
          }
          if (track.sign() < 0) {
            nPBar++;
          }
        }
        if (trackIsElectron) {
          if (idMethodEl == 0)
            fillIdentificationQA<recoAnalysisDir, kEl, tpcId, true>(recoAnalysis, track);
          if (idMethodEl == 1)
            fillIdentificationQA<recoAnalysisDir, kEl, tofId, true>(recoAnalysis, track);

          if (track.sign() > 0) {
            nElPlus++;
          }
          if (track.sign() < 0) {
            nElMinus++;
          }
        }
        if (trackIsDeuteron) {
          if (idMethodDe == 0)
            fillIdentificationQA<recoAnalysisDir, kDe, tpcId, true>(recoAnalysis, track);
          if (idMethodDe == 1)
            fillIdentificationQA<recoAnalysisDir, kDe, tofId, true>(recoAnalysis, track);

          if (track.sign() > 0) {
            nDePlus++;
          }
          if (track.sign() < 0) {
            nDeMinus++;
          }
        }
        recoAnalysis.fill(HIST("recoAnalysis/SelectedTrack_IdentificationTag"), trackIdTag);
      } // track loop ends

      nKaon = nKaPlus + nKaMinus;

      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_K0sPiKa"),
                        centFT0C, nTrack, nK0s,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nPiPlus, nPiMinus, nKaPlus, nKaMinus);
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_K0sPrDe"),
                        centFT0C, nTrack, nK0s,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nProton, nPBar, nDePlus, nDeMinus);
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_K0sKaEl"),
                        centFT0C, nTrack, nK0s, nRejectedPiPlus, nRejectedPiMinus,
                        nKaPlus, nKaMinus, nElPlus, nElMinus);
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_PiKaPr"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nPiPlus, nPiMinus, nKaPlus, nKaMinus, nProton, nPBar);
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_PiElDe"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nPiPlus, nPiMinus, nElPlus, nElMinus, nDePlus, nDeMinus);
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_KaPrDe"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nKaPlus, nKaMinus, nProton, nPBar, nDePlus, nDeMinus);
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_PrElDe"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nProton, nPBar, nElPlus, nElMinus, nDePlus, nDeMinus);

      recoAnalysis.fill(HIST("recoAnalysis/Sparse_newDynm_K0s_Ka"),
                        centFT0C, nTrack, nK0s, nKaon,
                        nK0s * nK0s, nKaon * nKaon, nK0s * nKaon);
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_newDynm_Kp_Km"),
                        centFT0C, nTrack, nKaPlus, nKaMinus,
                        nKaPlus * nKaPlus, nKaMinus * nKaMinus, nKaPlus * nKaMinus);

      // Collisions QA
      recoEvent.fill(HIST("recoEvent/h01_CollisionCount"), 0.5);
      recoEvent.fill(HIST("recoEvent/h02_VertexXRec"), collision.posX());
      recoEvent.fill(HIST("recoEvent/h03_VertexYRec"), collision.posY());
      recoEvent.fill(HIST("recoEvent/h04_VertexZRec"), collision.posZ());
      recoEvent.fill(HIST("recoEvent/h05_Centrality"), collision.centFT0C());
      recoEvent.fill(HIST("recoEvent/h06_V0Size"), v0sTablePerColl.size());
      recoEvent.fill(HIST("recoEvent/h07_TracksSize"), tracksTablePerColl.size());
      recoEvent.fill(HIST("recoEvent/h08_nTrack"), nTrack);
      recoEvent.fill(HIST("recoEvent/h09_nK0s"), nK0s);
      recoEvent.fill(HIST("recoEvent/h10_nPiPlus"), nPiPlus);
      recoEvent.fill(HIST("recoEvent/h11_nPiMinus"), nPiMinus);
      recoEvent.fill(HIST("recoEvent/h12_nKaPlus"), nKaPlus);
      recoEvent.fill(HIST("recoEvent/h13_nKaMinus"), nKaMinus);
      recoEvent.fill(HIST("recoEvent/h14_nProton"), nProton);
      recoEvent.fill(HIST("recoEvent/h15_nPBar"), nPBar);
      recoEvent.fill(HIST("recoEvent/h16_nElPlus"), nElPlus);
      recoEvent.fill(HIST("recoEvent/h17_nElMinus"), nElMinus);
      recoEvent.fill(HIST("recoEvent/h18_nDePlus"), nDePlus);
      recoEvent.fill(HIST("recoEvent/h19_nDeMinus"), nDeMinus);
    } // collision loop ends
  } // Process Function Ends
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<KaonIsospinFluctuations>(cfgc)};
}
