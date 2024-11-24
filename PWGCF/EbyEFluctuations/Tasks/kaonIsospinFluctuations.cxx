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
/// \brief  Kaon Isospin fluctuations
///
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Sadhana Dash (sadhana@phy.iitb.ac.in)

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

struct isospin_fluctuation {
  // Hisogram redistry:
  HistogramRegistry recoEvent{"recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoV0s{"recoV0s", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoK0s{"recoK0s", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoTracks{"recoTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoAnalysis{"recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoLambda{"recoLambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Configurables
  // Event Selection
  Configurable<float> cutZvertex{"cutZvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos to PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg to PV"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcaV0dau", 1, "DCA V0 Daughters"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};

  // Configurable K0s  //0.480-0.515 GeV/c2
  Configurable<double> mLowK0s{"mLowK0s", 0.48, "mLowK0s"};
  Configurable<double> mHighK0s{"mHighK0s", 0.515, "mHighK0s"};

  void init(InitContext const&)
  {
    // Print output histograms statistics
    // LOG(info) << "Starting Init";
    LOGF(info, "Starting init");

    std::vector<double> AxisSquarredArray;
    AxisSquarredArray.push_back(-1.5);
    AxisSquarredArray.push_back(-0.5);
    for (int i = 1; i < 1000; i++) {
      AxisSquarredArray.push_back(double(i * i) - 0.5);
      AxisSquarredArray.push_back(double(i * i) + 0.5);
    }

    // Axes
    // const AxisSpec axisCounts{1, 0., 1., ""};
    AxisSpec Axis_K0sMass = {200, 0.40f, 0.60f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec Axis_LambdaMass = {200, 1.f, 1.2F, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec Axis_XiMass = {200, 1.28f, 1.36f, "#it{M}_{inv} [GEv/#it{c}^{2}]"};

    AxisSpec Axis_vertexZ = {30, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec Axis_centFT0C = {1020, -1.0, 101.0, "centFT0C(percentile)"};

    AxisSpec Axis_p = {200, 0.0f, 10.0f, "#it{p} (GeV/#it{c})"};
    AxisSpec Axis_pt = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec Axis_tpcInnerParam = {200, 0.0f, 10.0f, "#it{p}_{tpcInnerParam} (GeV/#it{c})"};
    AxisSpec Axis_tofExpMom = {200, 0.0f, 10.0f, "#it{p}_{tofExpMom} (GeV/#it{c})"};

    AxisSpec Axis_eta = {100, -5, 5, "#eta"};
    AxisSpec Axis_phi = {90, -1, 8, "#phi (radians)"};
    AxisSpec Axis_rapidity = {200, -10, 10, "Rapidity (y)"};
    AxisSpec Axis_dcaXY = {2000, -100, 100, "dcaXY"};
    AxisSpec Axis_dcaZ = {2000, -100, 100, "dcaZ"};
    AxisSpec Axis_Sign = {10, -5, 5, "track.sign"};

    AxisSpec Axis_Mult = {1200, -10, 110};
    AxisSpec Axis_particleCount = {510, -10, 500};
    AxisSpec SquareCountAxis1 = {250000, -1, 249999, "SquaredCountAxis"};
    AxisSpec VectorAxis = AxisSquarredArray;

    AxisSpec Axis_tpcSignal = {10010, -1, 1000, "tpcSignal"};
    AxisSpec Axis_tofBeta = {400, -2.0, 2.0, "tofBeta"};

    AxisSpec Axis_tpcNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TPC}^{Pi}"};
    AxisSpec Axis_tofNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TOF}^{Pi}"};
    AxisSpec Axis_tpcNSigmaKa = {200, -10.0, 10.0, "n#sigma_{TPC}^{Ka}"};
    AxisSpec Axis_tofNSigmaKa = {200, -10.0, 10.0, "n#sigma_{TOF}^{Ka}"};
    AxisSpec Axis_tpcNSigmaPr = {200, -10.0, 10.0, "n#sigma_{TPC}^{Pr}"};
    AxisSpec Axis_tofNSigmaPr = {200, -10.0, 10.0, "n#sigma_{TOF}^{Pr}"};
    AxisSpec Axis_tpcNSigmaEl = {200, -10.0, 10.0, "n#sigma_{TPC}^{El}"};
    AxisSpec Axis_tofNSigmaEl = {200, -10.0, 10.0, "n#sigma_{TOF}^{El}"};
    AxisSpec Axis_tpcNSigmaDe = {200, -10.0, 10.0, "n#sigma_{TPC}^{De}"};
    AxisSpec Axis_tofNSigmaDe = {200, -10.0, 10.0, "n#sigma_{TOF}^{De}"};

    AxisSpec Axis_tpcNClsCrossedRows = {200, -1.5, 198.5, "tpcNClsCrossedRows"};
    AxisSpec Axis_isGlobalTrack = {4, -1, 3, "isGobalTrack"};
    AxisSpec Axis_isK0sDau = {4, -1, 3, "isK0sDau"};

    AxisSpec Axis_00 = Axis_p;             // p              //0
    AxisSpec Axis_01 = Axis_pt;            // pt             //1
    AxisSpec Axis_02 = Axis_tpcInnerParam; // tpcInnerParam  //2
    AxisSpec Axis_03 = Axis_tofExpMom;     // tofExpMom      //3

    AxisSpec Axis_05 = Axis_tpcSignal; // tpcSignal      //5
    AxisSpec Axis_06 = Axis_tofBeta;   // Axis_tofBeta   //

    AxisSpec Axis_20 = Axis_tpcNSigmaPi; // tpcNSigmaPi    //5
    AxisSpec Axis_21 = Axis_tofNSigmaPi; // tofNSigmaPi
    AxisSpec Axis_22 = Axis_tpcNSigmaKa; // tpcNSigmaKa    //5
    AxisSpec Axis_23 = Axis_tofNSigmaKa; // tofNSigmaKa
    AxisSpec Axis_24 = Axis_tpcNSigmaPr; // tpcNSigmaPr    //5
    AxisSpec Axis_25 = Axis_tofNSigmaPr; // tofNSigmaPr
    AxisSpec Axis_26 = Axis_tpcNSigmaEl; // tpcNSigmaEl    //5
    AxisSpec Axis_27 = Axis_tofNSigmaEl; // tofNSigmaEl
    AxisSpec Axis_28 = Axis_tpcNSigmaDe; // tpcNSigmaDe    //5
    AxisSpec Axis_29 = Axis_tofNSigmaDe; // tofNSigmaDe

    // Histograms
    // V0 Information
    recoV0s.add("hV0s_0_0_00_hV0TagCount", "hV0s_0_0_00_hV0TagCount", {HistType::kTH1F, {{12, -2, 10}}});         // 001 = Kaon, 010 = Lambda, 100 = AnitLambda
    recoV0s.add("hV0s_0_0_00_hTrueV0TagCount", "hV0s_0_0_00_hTrueV0TagCount", {HistType::kTH1F, {{12, -2, 10}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda

    recoV0s.add("hV0s_0_0_01_K0s_Mass_Full", "hV0s_0_0_01_K0s_Mass_Full", {HistType::kTH1F, {Axis_K0sMass}});
    recoV0s.add("hV0s_0_0_02_Lambda_Mass_Full", "hV0s_0_0_02_Lambda_Mass_Full", {HistType::kTH1F, {Axis_LambdaMass}});
    recoV0s.add("hV0s_0_0_03_AntiLambda_Mass_Full", "hV0s_0_0_03_AntiLambda_Mass_Full", {HistType::kTH1F, {Axis_LambdaMass}});

    // K0s topological/PID cuts
    recoV0s.add("hV0s_0_1_01_K0s_dcapostopv", "hV0s_K0s_dcapostopv", {HistType::kTH1F, {{100000, -50, 50}}});
    recoV0s.add("hV0s_0_1_02_K0s_dcanegtopv", "hV0s_K0s_dcanegtopv", {HistType::kTH1F, {{100000, -50, 50}}});
    recoV0s.add("hV0s_0_1_03_K0s_dcaV0daughters", "hV0s_K0s_dcaV0daughters", {HistType::kTH1F, {{2000, -10.0, 10.0}}});
    recoV0s.add("hV0s_0_1_04_K0s_v0cosPA", "hV0s_K0s_v0cosPA", {HistType::kTH1F, {{3000, -1.5, 1.5}}});
    recoV0s.add("hV0s_0_1_05_K0s_v0radius", "hV0s_K0s_v0radius", {HistType::kTH1F, {{100000, -50, 50}}});

    // K0s-FullInformation
    recoV0s.add("hV0s_table_0_2_01_K0s_Mass_Selected", "hV0s_table_0_2_01_K0s_Mass_Selected", {HistType::kTH1F, {Axis_K0sMass}});
    recoV0s.add("hV0s_table_0_2_02_K0s_P", "hV0s_table_0_2_02_K0s_P", {HistType::kTH1F, {Axis_p}});
    recoV0s.add("hV0s_table_0_2_03_K0s_Pt", "hV0s_table_0_2_03_K0s_Pt", {HistType::kTH1F, {Axis_pt}});
    recoV0s.add("hV0s_table_0_2_04_K0s_Eta", "hV0s_table_0_2_04_K0s_Eta", {HistType::kTH1F, {Axis_eta}});
    recoV0s.add("hV0s_table_0_2_05_K0s_Phi", "hV0s_table_0_2_05_K0s_Phi", {HistType::kTH1F, {Axis_phi}});
    recoV0s.add("hV0s_table_0_2_06_K0s_Rapidity", "hV0s_table_0_2_06_K0s_Rapidity", {HistType::kTH1F, {Axis_rapidity}});

    // Pion identification from K0s
    recoV0s.add("hV0s_table_0_3_01_K0s_Pi_P", "hV0s_table_0_3_01_K0s_Pi_P", kTH1F, {Axis_p});
    recoV0s.add("hV0s_table_0_3_01_K0s_Pi_Pt", "hV0s_table_0_3_01_K0s_Pi_Pt", kTH1F, {Axis_pt});                                  // Axis_pt            ;
    recoV0s.add("hV0s_table_0_3_01_K0s_Pi_tpcInnerParam", "hV0s_table_0_3_01_K0s_Pi_tpcInnerParam", kTH1F, {Axis_tpcInnerParam}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_01_K0s_Pi_tofExpMom", "hV0s_table_0_3_01_K0s_Pi_tofExpMom", kTH1F, {Axis_tofExpMom});             // Axis_tofExpMom     ;

    recoV0s.add("hV0s_table_0_3_02_K0s_Pi_Eta", "hV0s_table_0_3_02_K0s_Pi_Eta", {HistType::kTH1F, {Axis_eta}});
    recoV0s.add("hV0s_table_0_3_03_K0s_Pi_Phi", "hV0s_table_0_3_03_K0s_Pi_Phi", {HistType::kTH1F, {Axis_phi}});
    recoV0s.add("hV0s_table_0_3_04_K0s_Pi_Rapidity", "hV0s_table_0_3_04_K0s_Pi_Rapidity", {HistType::kTH1F, {Axis_rapidity}});
    recoV0s.add("hV0s_table_0_3_05_K0s_Pi_isPVContributor", "hV0s_table_0_3_05_K0s_Pi_isPVContributor", kTH1D, {{4, -1, 3}});
    recoV0s.add("hV0s_table_0_3_06_K0s_Pi_isGlobalTrack", "hV0s_table_0_3_06_K0s_Pi_isGlobalTrack", kTH1D, {{4, -1, 3}});
    recoV0s.add("hV0s_table_0_3_07_K0s_Pi_DcaXY", "hV0s_table_0_3_07_K0s_Pi_DcaXY", kTH1D, {Axis_dcaXY});
    recoV0s.add("hV0s_table_0_3_08_K0s_Pi_DcaZ", "hV0s_table_0_3_08_K0s_Pi_DcaZ", kTH1D, {Axis_dcaZ});

    // momemtum
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_01", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_02", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_03", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;
    // tpcSignal
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_05", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_05", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_05", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;
    // tofBeta
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_06", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_06", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_06", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;
    // Look at Pion
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_20", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_20", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_20", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_20", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_21", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_21", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_21", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_21", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_20_21", "hV0s_table_0_3_10_K0s_Pi_Id0_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;

    // momemtum
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_01", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_02", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_03", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;
    // tpcSignal
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_05", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_05", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_05", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;
    // tofBeta
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_06", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_06", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_06", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;
    // Look at Pion
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_20", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_20", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_20", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_20", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_21", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_21", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_21", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_21", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_20_21", "hV0s_table_0_3_10_K0s_Pi_Id1_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;
    // K0s-FullInformation

    // K0s-Selected
    //  Mass
    recoV0s.add("hV0s_0_2_01_K0s_Mass_DynamicCuts", "hV0s_K0s_Mass_DynamicCuts", {HistType::kTH1F, {Axis_K0sMass}});

    // Pt, eta, phi, rapidity of K0s
    recoV0s.add("hV0s_0_2_01_K0s_Mass_Selected_LambdaMass", "hV0s_0_2_01_K0s_Mass_Selected_LambdaMass", {HistType::kTH1F, {Axis_LambdaMass}});
    recoV0s.add("hV0s_0_2_01_K0s_Mass_Selected_AntiLambdaMass", "hV0s_0_2_01_K0s_Mass_Selected_AntiLambdaMass", {HistType::kTH1F, {Axis_LambdaMass}});
    recoV0s.add("hV0s_0_2_01_K0s_Mass_Selected", "hV0s_0_2_01_K0s_Mass_Selected", {HistType::kTH1F, {Axis_K0sMass}});
    recoV0s.add("hV0s_0_2_02_K0s_P", "hV0s_0_2_02_K0s_P", {HistType::kTH1F, {Axis_p}});
    recoV0s.add("hV0s_0_2_03_K0s_Pt", "hV0s_0_2_03_K0s_Pt", {HistType::kTH1F, {Axis_pt}});
    recoV0s.add("hV0s_0_2_04_K0s_Eta", "hV0s_0_2_04_K0s_Eta", {HistType::kTH1F, {Axis_eta}});
    recoV0s.add("hV0s_0_2_05_K0s_Phi", "hV0s_0_2_05_K0s_Phi", {HistType::kTH1F, {Axis_phi}});
    recoV0s.add("hV0s_0_2_06_K0s_Rapidity", "hV0s_0_2_06_K0s_Rapidity", {HistType::kTH1F, {Axis_rapidity}});

    // Pion identification from K0s
    recoV0s.add("hV0s_0_3_01_K0s_Pi_P", "hV0s_0_3_01_K0s_Pi_P", kTH1F, {Axis_p});
    recoV0s.add("hV0s_0_3_01_K0s_Pi_Pt", "hV0s_0_3_01_K0s_Pi_Pt", kTH1F, {Axis_pt});                                  // Axis_pt            ;
    recoV0s.add("hV0s_0_3_01_K0s_Pi_tpcInnerParam", "hV0s_0_3_01_K0s_Pi_tpcInnerParam", kTH1F, {Axis_tpcInnerParam}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_01_K0s_Pi_tofExpMom", "hV0s_0_3_01_K0s_Pi_tofExpMom", kTH1F, {Axis_tofExpMom});             // Axis_tofExpMom     ;

    recoV0s.add("hV0s_0_3_02_K0s_Pi_Eta", "hV0s_0_3_02_K0s_Pi_Eta", {HistType::kTH1F, {Axis_eta}});
    recoV0s.add("hV0s_0_3_03_K0s_Pi_Phi", "hV0s_0_3_03_K0s_Pi_Phi", {HistType::kTH1F, {Axis_phi}});
    recoV0s.add("hV0s_0_3_04_K0s_Pi_Rapidity", "hV0s_0_3_04_K0s_Pi_Rapidity", {HistType::kTH1F, {Axis_rapidity}});
    recoV0s.add("hV0s_0_3_05_K0s_Pi_isPVContributor", "hV0s_0_3_05_K0s_Pi_isPVContributor", kTH1D, {{4, -1, 3}});
    recoV0s.add("hV0s_0_3_06_K0s_Pi_isGlobalTrack", "hV0s_0_3_06_K0s_Pi_isGlobalTrack", kTH1D, {{4, -1, 3}});
    recoV0s.add("hV0s_0_3_07_K0s_Pi_DcaXY", "hV0s_0_3_07_K0s_Pi_DcaXY", kTH1D, {Axis_dcaXY});
    recoV0s.add("hV0s_0_3_08_K0s_Pi_DcaZ", "hV0s_0_3_08_K0s_Pi_DcaZ", kTH1D, {Axis_dcaZ});

    // momemtum
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_01", "hV0s_0_3_10_K0s_Pi_Id0_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_02", "hV0s_0_3_10_K0s_Pi_Id0_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_03", "hV0s_0_3_10_K0s_Pi_Id0_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;
    // tpcSignal
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_05", "hV0s_0_3_10_K0s_Pi_Id0_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_05", "hV0s_0_3_10_K0s_Pi_Id0_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_05", "hV0s_0_3_10_K0s_Pi_Id0_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;
    // tofBeta
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_06", "hV0s_0_3_10_K0s_Pi_Id0_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_06", "hV0s_0_3_10_K0s_Pi_Id0_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_06", "hV0s_0_3_10_K0s_Pi_Id0_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;
    // Look at Pion
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_20", "hV0s_0_3_10_K0s_Pi_Id0_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_01_20", "hV0s_0_3_10_K0s_Pi_Id0_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_20", "hV0s_0_3_10_K0s_Pi_Id0_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_20", "hV0s_0_3_10_K0s_Pi_Id0_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_21", "hV0s_0_3_10_K0s_Pi_Id0_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_01_21", "hV0s_0_3_10_K0s_Pi_Id0_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_21", "hV0s_0_3_10_K0s_Pi_Id0_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_21", "hV0s_0_3_10_K0s_Pi_Id0_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id0_Axis_20_21", "hV0s_0_3_10_K0s_Pi_Id0_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;

    // momemtum
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_01", "hV0s_0_3_10_K0s_Pi_Id1_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_02", "hV0s_0_3_10_K0s_Pi_Id1_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_03", "hV0s_0_3_10_K0s_Pi_Id1_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;
    // tpcSignal
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_05", "hV0s_0_3_10_K0s_Pi_Id1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_05", "hV0s_0_3_10_K0s_Pi_Id1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_05", "hV0s_0_3_10_K0s_Pi_Id1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;
    // tofBeta
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_06", "hV0s_0_3_10_K0s_Pi_Id1_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_06", "hV0s_0_3_10_K0s_Pi_Id1_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_06", "hV0s_0_3_10_K0s_Pi_Id1_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;
    // Look at Pion
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_20", "hV0s_0_3_10_K0s_Pi_Id1_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_01_20", "hV0s_0_3_10_K0s_Pi_Id1_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_20", "hV0s_0_3_10_K0s_Pi_Id1_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_20", "hV0s_0_3_10_K0s_Pi_Id1_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_21", "hV0s_0_3_10_K0s_Pi_Id1_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_01_21", "hV0s_0_3_10_K0s_Pi_Id1_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_21", "hV0s_0_3_10_K0s_Pi_Id1_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_21", "hV0s_0_3_10_K0s_Pi_Id1_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoV0s.add("hV0s_0_3_10_K0s_Pi_Id1_Axis_20_21", "hV0s_0_3_10_K0s_Pi_Id1_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;

    recoV0s.add("hV0s_0_4_01_K0s_v0DaughterCollisionIndexTag", "hV0s_K0s_v0DaughterCollisionIndexTag", {HistType::kTH1D, {{22, -1.0, 10.0}}});
    recoV0s.add("hV0s_0_4_02_K0s_nCommonPionOfDifferentK0s", "hV0s_K0s_nCommonPionOfDifferentK0s", {HistType::kTH1D, {{44, -2, 20}}});
    // K0s-Selected
    //
    // Event Selection
    recoEvent.add("hCollisionCount", "hCollisionCount", {HistType::kTH1D, {{1, 0, 1}}});
    recoEvent.add("hVertexXRec", "hVertexXRec", {HistType::kTH1D, {{10000, -0.2, 0.2}}});
    recoEvent.add("hVertexYRec", "hVertexYRec", {HistType::kTH1D, {{10000, -0.2, 0.2}}});
    recoEvent.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {Axis_vertexZ}});
    recoEvent.add("hV0Size", "hV0Size", {HistType::kTH1F, {{60, -10, 50}}});
    recoEvent.add("hCascadeSize", "hCascadeSize", {HistType::kTH1F, {{60, -10, 50}}});
    recoEvent.add("hCentrality", "hCentrality", {HistType::kTH1F, {Axis_centFT0C}});
    //
    // K0s reconstruction
    // K0s topological/PID cuts
    recoK0s.add("hK0s_1_01_dcapostopv", "hK0s_dcapostopv", {HistType::kTH1F, {{100000, -50, 50}}});
    recoK0s.add("hK0s_1_02_dcanegtopv", "hK0s_dcanegtopv", {HistType::kTH1F, {{100000, -50, 50}}});
    recoK0s.add("hK0s_1_03_dcaV0daughters", "hK0s_dcaV0daughters", {HistType::kTH1F, {{2000, -10.0, 10.0}}});
    recoK0s.add("hK0s_1_04_v0cosPA", "hK0s_v0cosPA", {HistType::kTH1F, {{3000, -1.5, 1.5}}});
    recoK0s.add("hK0s_1_05_v0radius", "hK0s_v0radius", {HistType::kTH1F, {{100000, -50, 50}}});
    // Mass
    recoK0s.add("hK0s_2_00_Mass_Full", "hK0s_Mass_Full", {HistType::kTH1F, {Axis_K0sMass}});
    recoK0s.add("hK0s_2_01_Mass_DynamicCuts", "hK0s_Mass_DynamicCuts", {HistType::kTH1F, {Axis_K0sMass}});

    // Pt, eta, phi, rapidity of K0s
    recoK0s.add("hK0s_2_01_K0s_Mass_Selected_LambdaMass", "hK0s_2_01_K0s_Mass_Selected_LambdaMass", {HistType::kTH1F, {Axis_LambdaMass}});
    recoK0s.add("hK0s_2_01_K0s_Mass_Selected_AntiLambdaMass", "hK0s_2_01_K0s_Mass_Selected_AntiLambdaMass", {HistType::kTH1F, {Axis_LambdaMass}});
    recoK0s.add("hK0s_2_01_Mass_Selected", "hK0s_Mass_Selected", {HistType::kTH1F, {Axis_K0sMass}});
    recoK0s.add("hK0s_2_02_P", "hK0s_2_02_P", {HistType::kTH1F, {Axis_p}});
    recoK0s.add("hK0s_2_03_Pt", "hK0s_Pt", {HistType::kTH1F, {Axis_pt}});
    recoK0s.add("hK0s_2_04_Eta", "hK0s_Eta", {HistType::kTH1F, {Axis_eta}});
    recoK0s.add("hK0s_2_05_Phi", "hK0s_Phi", {HistType::kTH1F, {Axis_phi}});
    recoK0s.add("hK0s_2_06_Rapidity", "hK0s_Rapidity", {HistType::kTH1F, {Axis_rapidity}});

    // Pion identification from K0s
    recoK0s.add("hK0s_3_01_Pion_P", "hK0s_3_01_Pion_P", kTH1F, {Axis_p});
    recoK0s.add("hK0s_3_01_Pion_Pt", "hK0s_3_01_Pion_Pt", kTH1F, {Axis_pt});                                  // Axis_pt            ;
    recoK0s.add("hK0s_3_01_Pion_tpcInnerParam", "hK0s_3_01_Pion_tpcInnerParam", kTH1F, {Axis_tpcInnerParam}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_01_Pion_tofExpMom", "hK0s_3_01_Pion_tofExpMom", kTH1F, {Axis_tofExpMom});             // Axis_tofExpMom     ;

    recoK0s.add("hK0s_3_02_Pion_Eta", "hK0s_Pi_Eta", {HistType::kTH1F, {Axis_eta}});
    recoK0s.add("hK0s_3_03_Pion_Phi", "hK0s_Pi_Phi", {HistType::kTH1F, {Axis_phi}});
    recoK0s.add("hK0s_3_04_Pion_Rapidity", "hK0s_Pi_Rapidity", {HistType::kTH1F, {Axis_rapidity}});
    recoK0s.add("hK0s_3_05_Daughter_isPVContributor", "hK0s_Daughter_isPVContributor", kTH1D, {{4, -1, 3}});
    recoK0s.add("hK0s_3_06_Daughter_isGlobalTrack", "hK0s_Daughter_isGlobalTrack", kTH1D, {{4, -1, 3}});
    recoK0s.add("hK0s_3_07_Daughter_DcaXY", "hK0s_Daughter_DcaXY", kTH1D, {Axis_dcaXY});
    recoK0s.add("hK0s_3_08_Daughter_DcaZ", "hK0s_Daughter_DcaZ", kTH1D, {Axis_dcaZ});

    // momemtum
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_00_01", "htrack_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_00_02", "htrack_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_00_03", "htrack_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;
    // tpcSignal
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_00_05", "htrack_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_02_05", "htrack_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_03_05", "htrack_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;
    // tofBeta
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_00_06", "htrack_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_02_06", "htrack_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_03_06", "htrack_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;
    // Look at Pion
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_00_20", "htrack_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_01_20", "htrack_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_02_20", "htrack_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_03_20", "htrack_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_00_21", "htrack_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_01_21", "htrack_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_02_21", "htrack_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_03_21", "htrack_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id0_Axis_20_21", "htrack_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;

    // momemtum
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_00_01", "htrack_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_00_02", "htrack_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_00_03", "htrack_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;
    // tpcSignal
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_00_05", "htrack_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_02_05", "htrack_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_03_05", "htrack_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;
    // tofBeta
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_00_06", "htrack_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_02_06", "htrack_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_03_06", "htrack_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;
    // Look at Pion
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_00_20", "htrack_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_01_20", "htrack_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_02_20", "htrack_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_03_20", "htrack_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_00_21", "htrack_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_01_21", "htrack_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_02_21", "htrack_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_03_21", "htrack_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoK0s.add("hK0s_3_10_K0s_Pi_Id1_Axis_20_21", "htrack_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;

    recoK0s.add("hK0s_4_01_v0DaughterCollisionIndexTag", "hK0s_v0DaughterCollisionIndexTag", {HistType::kTH1D, {{9, -1, 10}}});

    recoK0s.add("hK0s_5_01_centFTOC_mK0s", "hK0s_5_01_centFTOC_mK0s", kTH2F, {Axis_centFT0C, Axis_K0sMass});
    //
    // Lambda reconstruction
    // // Mass
    recoLambda.add("hLambda_Mass_Full", "hLambda_Mass_Full", {HistType::kTH1F, {Axis_LambdaMass}});
    recoLambda.add("hAntiLambda_Mass_Full", "hAntiLambda_Mass_Full", {HistType::kTH1F, {Axis_LambdaMass}});
    // Tracks reconstruction
    // FullTrack
    recoTracks.add("htracks_00_1_0_FullTrack_P", "hTracks_SelectedTrack_P", {HistType::kTH1F, {Axis_p}});
    recoTracks.add("htracks_00_1_0_FullTrack_tpcInnerParam", "hTracks_SelectedTrack_tpcInnerParam", {HistType::kTH1F, {Axis_tpcInnerParam}});
    recoTracks.add("htracks_00_1_0_FullTrack_tofExpMom", "hTracks_SelectedTrack_tofExpMom", {HistType::kTH1F, {Axis_tofExpMom}});

    recoTracks.add("hTracks_00_1_1_FullTrack_Pt", "hTracks_FullTrack_Pt", {HistType::kTH1F, {Axis_pt}});
    recoTracks.add("hTracks_00_1_2_FullTrack_Eta", "hTracks_FullTrack_Eta", {HistType::kTH1F, {Axis_eta}});
    recoTracks.add("hTracks_00_1_3_FullTrack_Phi", "hTracks_FullTrack_Phi", {HistType::kTH1F, {Axis_phi}});
    recoTracks.add("hTracks_00_1_4_FullTrack_DcaXY", "hTracks_FullTrack_DcaXY", {HistType::kTH1F, {Axis_dcaXY}});
    recoTracks.add("hTracks_00_1_5_FullTrack_DcaZ", "hTracks_FullTrack_DcaZ", {HistType::kTH1F, {Axis_dcaZ}});
    recoTracks.add("hTracks_00_1_6_FullTrack_Sign", "hTracks_FullTrack_Sign", {HistType::kTH1D, {Axis_Sign}});

    // DcaXY
    recoTracks.add("htracks_00_1_7_FullTrack_00_DcaXY", "htracks_FullTrack_00_DcaXY", kTH2F, {Axis_00, Axis_dcaXY});
    recoTracks.add("htracks_00_1_7_FullTrack_01_DcaXY", "htracks_FullTrack_01_DcaXY", kTH2F, {Axis_01, Axis_dcaXY});
    recoTracks.add("htracks_00_1_7_FullTrack_02_DcaXY", "htracks_FullTrack_02_DcaXY", kTH2F, {Axis_02, Axis_dcaXY});
    recoTracks.add("htracks_00_1_7_FullTrack_03_DcaXY", "htracks_FullTrack_03_DcaXY", kTH2F, {Axis_03, Axis_dcaXY});

    // DcaZ
    recoTracks.add("htracks_00_1_7_FullTrack_00_DcaZ", "htracks_FullTrack_00_DcaZ", kTH2F, {Axis_00, Axis_dcaZ});
    recoTracks.add("htracks_00_1_7_FullTrack_01_DcaZ", "htracks_FullTrack_01_DcaZ", kTH2F, {Axis_01, Axis_dcaZ});
    recoTracks.add("htracks_00_1_7_FullTrack_02_DcaZ", "htracks_FullTrack_02_DcaZ", kTH2F, {Axis_02, Axis_dcaZ});
    recoTracks.add("htracks_00_1_7_FullTrack_03_DcaZ", "htracks_FullTrack_03_DcaZ", kTH2F, {Axis_03, Axis_dcaZ});

    // momemtum
    recoTracks.add("hTracks_00_2_1_FullTrack_Axis_00_01", "htrack_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoTracks.add("hTracks_00_2_2_FullTrack_Axis_00_02", "htrack_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_2_3_FullTrack_Axis_00_03", "htrack_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;

    // tpcSignal
    recoTracks.add("hTracks_00_3_1_FullTrack_Axis_00_05", "htrack_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoTracks.add("hTracks_00_3_2_FullTrack_Axis_02_05", "htrack_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_3_3_FullTrack_Axis_03_05", "htrack_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;

    // tofBeta
    recoTracks.add("hTracks_00_4_1_FullTrack_Axis_00_06", "htrack_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoTracks.add("hTracks_00_4_2_FullTrack_Axis_02_06", "htrack_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_4_3_FullTrack_Axis_03_06", "htrack_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;

    // Look at Pion
    recoTracks.add("hTracks_00_5_1_FullTrack_Axis_00_20", "htrack_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoTracks.add("hTracks_00_5_2_FullTrack_Axis_01_20", "htrack_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoTracks.add("hTracks_00_5_3_FullTrack_Axis_02_20", "htrack_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_5_4_FullTrack_Axis_03_20", "htrack_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_5_5_FullTrack_Axis_00_21", "htrack_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoTracks.add("hTracks_00_5_6_FullTrack_Axis_01_21", "htrack_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoTracks.add("hTracks_00_5_7_FullTrack_Axis_02_21", "htrack_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_5_8_FullTrack_Axis_03_21", "htrack_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_5_9_FullTrack_Axis_20_21", "htrack_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    recoTracks.add("hTracks_00_6_1_FullTrack_Axis_00_22", "htrack_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // Axis_p             ;
    recoTracks.add("hTracks_00_6_2_FullTrack_Axis_01_22", "htrack_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // Axis_pt            ;
    recoTracks.add("hTracks_00_6_3_FullTrack_Axis_02_22", "htrack_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_6_4_FullTrack_Axis_03_22", "htrack_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_6_5_FullTrack_Axis_00_23", "htrack_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // Axis_p             ;
    recoTracks.add("hTracks_00_6_6_FullTrack_Axis_01_23", "htrack_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // Axis_pt            ;
    recoTracks.add("hTracks_00_6_7_FullTrack_Axis_02_23", "htrack_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_6_8_FullTrack_Axis_03_23", "htrack_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_6_9_FullTrack_Axis_22_23", "htrack_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    recoTracks.add("hTracks_00_7_1_FullTrack_Axis_00_24", "htrack_Axis_00_24", kTH2F, {Axis_00, Axis_24}); // Axis_p             ;
    recoTracks.add("hTracks_00_7_2_FullTrack_Axis_01_24", "htrack_Axis_01_24", kTH2F, {Axis_01, Axis_24}); // Axis_pt            ;
    recoTracks.add("hTracks_00_7_3_FullTrack_Axis_02_24", "htrack_Axis_02_24", kTH2F, {Axis_02, Axis_24}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_7_4_FullTrack_Axis_03_24", "htrack_Axis_03_24", kTH2F, {Axis_03, Axis_24}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_7_5_FullTrack_Axis_00_25", "htrack_Axis_00_25", kTH2F, {Axis_00, Axis_25}); // Axis_p             ;
    recoTracks.add("hTracks_00_7_6_FullTrack_Axis_01_25", "htrack_Axis_01_25", kTH2F, {Axis_01, Axis_25}); // Axis_pt            ;
    recoTracks.add("hTracks_00_7_7_FullTrack_Axis_02_25", "htrack_Axis_02_25", kTH2F, {Axis_02, Axis_25}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_7_8_FullTrack_Axis_03_25", "htrack_Axis_03_25", kTH2F, {Axis_03, Axis_25}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_7_9_FullTrack_Axis_24_25", "htrack_Axis_24_25", kTH2F, {Axis_24, Axis_25}); // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    recoTracks.add("hTracks_00_8_1_FullTrack_Axis_00_26", "htrack_Axis_00_26", kTH2F, {Axis_00, Axis_26}); // Axis_p             ;
    recoTracks.add("hTracks_00_8_2_FullTrack_Axis_01_26", "htrack_Axis_01_26", kTH2F, {Axis_01, Axis_26}); // Axis_pt            ;
    recoTracks.add("hTracks_00_8_3_FullTrack_Axis_02_26", "htrack_Axis_02_26", kTH2F, {Axis_02, Axis_26}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_8_4_FullTrack_Axis_03_26", "htrack_Axis_03_26", kTH2F, {Axis_03, Axis_26}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_8_5_FullTrack_Axis_00_27", "htrack_Axis_00_27", kTH2F, {Axis_00, Axis_27}); // Axis_p             ;
    recoTracks.add("hTracks_00_8_6_FullTrack_Axis_01_27", "htrack_Axis_01_27", kTH2F, {Axis_01, Axis_27}); // Axis_pt            ;
    recoTracks.add("hTracks_00_8_7_FullTrack_Axis_02_27", "htrack_Axis_02_27", kTH2F, {Axis_02, Axis_27}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_8_8_FullTrack_Axis_03_27", "htrack_Axis_03_27", kTH2F, {Axis_03, Axis_27}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_8_9_FullTrack_Axis_26_27", "htrack_Axis_26_27", kTH2F, {Axis_26, Axis_27}); // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    recoTracks.add("hTracks_00_9_1_FullTrack_Axis_00_28", "htrack_Axis_00_28", kTH2F, {Axis_00, Axis_28}); // Axis_p             ;
    recoTracks.add("hTracks_00_9_2_FullTrack_Axis_01_28", "htrack_Axis_01_28", kTH2F, {Axis_01, Axis_28}); // Axis_pt            ;
    recoTracks.add("hTracks_00_9_3_FullTrack_Axis_02_28", "htrack_Axis_02_28", kTH2F, {Axis_02, Axis_28}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_9_4_FullTrack_Axis_03_28", "htrack_Axis_03_28", kTH2F, {Axis_03, Axis_28}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_9_5_FullTrack_Axis_00_29", "htrack_Axis_00_29", kTH2F, {Axis_00, Axis_29}); // Axis_p             ;
    recoTracks.add("hTracks_00_9_6_FullTrack_Axis_01_29", "htrack_Axis_01_29", kTH2F, {Axis_01, Axis_29}); // Axis_pt            ;
    recoTracks.add("hTracks_00_9_7_FullTrack_Axis_02_29", "htrack_Axis_02_29", kTH2F, {Axis_02, Axis_29}); // Axis_tpcInnerParam ;
    recoTracks.add("hTracks_00_9_8_FullTrack_Axis_03_29", "htrack_Axis_03_29", kTH2F, {Axis_03, Axis_29}); // Axis_tofExpMom     ;
    recoTracks.add("hTracks_00_9_9_FullTrack_Axis_28_29", "htrack_Axis_28_29", kTH2F, {Axis_28, Axis_29}); // Axis_tpcInnerParam ;
    // Deuteron
    // FullTrack

    // SelectedTrack
    recoTracks.add("htracks_11_1_0_SelectedTrack_P", "hTracks_SelectedTrack_P", {HistType::kTH1F, {Axis_p}});
    recoTracks.add("htracks_11_1_0_SelectedTrack_tpcInnerParam", "hTracks_SelectedTrack_tpcInnerParam", {HistType::kTH1F, {Axis_tpcInnerParam}});
    recoTracks.add("htracks_11_1_0_SelectedTrack_tofExpMom", "hTracks_SelectedTrack_tofExpMom", {HistType::kTH1F, {Axis_tofExpMom}});

    recoTracks.add("htracks_11_1_1_SelectedTrack_Pt", "hTracks_SelectedTrack_Pt", {HistType::kTH1F, {Axis_pt}});
    recoTracks.add("htracks_11_1_2_SelectedTrack_Eta", "hTracks_SelectedTrack_Eta", {HistType::kTH1F, {Axis_eta}});
    recoTracks.add("htracks_11_1_3_SelectedTrack_Phi", "hTracks_SelectedTrack_Phi", {HistType::kTH1F, {Axis_phi}});
    recoTracks.add("htracks_11_1_4_SelectedTrack_DcaXY", "hTracks_SelectedTrack_DcaXY", {HistType::kTH1F, {Axis_dcaXY}});
    recoTracks.add("htracks_11_1_5_SelectedTrack_DcaZ", "hTracks_SelectedTrack_DcaZ", {HistType::kTH1F, {Axis_dcaZ}});
    recoTracks.add("htracks_11_1_6_SelectedTrack_Sign", "hTracks_SelectedTrack_Sign", {HistType::kTH1D, {Axis_Sign}});

    // DcaXY
    recoTracks.add("htracks_11_1_7_SelectedTrack_00_DcaXY", "htracks_SelectedTrack_00_DcaXY", kTH2F, {Axis_00, Axis_dcaXY});
    recoTracks.add("htracks_11_1_7_SelectedTrack_01_DcaXY", "htracks_SelectedTrack_01_DcaXY", kTH2F, {Axis_01, Axis_dcaXY});
    recoTracks.add("htracks_11_1_7_SelectedTrack_02_DcaXY", "htracks_SelectedTrack_02_DcaXY", kTH2F, {Axis_02, Axis_dcaXY});
    recoTracks.add("htracks_11_1_7_SelectedTrack_03_DcaXY", "htracks_SelectedTrack_03_DcaXY", kTH2F, {Axis_03, Axis_dcaXY});

    // DcaZ
    recoTracks.add("htracks_11_1_7_SelectedTrack_00_DcaZ", "htracks_SelectedTrack_00_DcaZ", kTH2F, {Axis_00, Axis_dcaZ});
    recoTracks.add("htracks_11_1_7_SelectedTrack_01_DcaZ", "htracks_SelectedTrack_01_DcaZ", kTH2F, {Axis_01, Axis_dcaZ});
    recoTracks.add("htracks_11_1_7_SelectedTrack_02_DcaZ", "htracks_SelectedTrack_02_DcaZ", kTH2F, {Axis_02, Axis_dcaZ});
    recoTracks.add("htracks_11_1_7_SelectedTrack_03_DcaZ", "htracks_SelectedTrack_03_DcaZ", kTH2F, {Axis_03, Axis_dcaZ});

    // momemtum
    recoTracks.add("htracks_11_2_1_SelectedTrack_Axis_00_01", "htrack_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    recoTracks.add("htracks_11_2_2_SelectedTrack_Axis_00_02", "htrack_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_2_3_SelectedTrack_Axis_00_03", "htrack_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;

    // tpcSignal
    recoTracks.add("htracks_11_3_1_SelectedTrack_Axis_00_05", "htrack_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    recoTracks.add("htracks_11_3_2_SelectedTrack_Axis_02_05", "htrack_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_3_3_SelectedTrack_Axis_03_05", "htrack_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;

    // tofBeta
    recoTracks.add("htracks_11_4_1_SelectedTrack_Axis_00_06", "htrack_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    recoTracks.add("htracks_11_4_2_SelectedTrack_Axis_02_06", "htrack_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_4_3_SelectedTrack_Axis_03_06", "htrack_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;

    // Look at Pion
    recoTracks.add("htracks_11_5_1_SelectedTrack_Axis_00_20", "htrack_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    recoTracks.add("htracks_11_5_2_SelectedTrack_Axis_01_20", "htrack_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    recoTracks.add("htracks_11_5_3_SelectedTrack_Axis_02_20", "htrack_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_5_4_SelectedTrack_Axis_03_20", "htrack_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_5_5_SelectedTrack_Axis_00_21", "htrack_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    recoTracks.add("htracks_11_5_6_SelectedTrack_Axis_01_21", "htrack_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    recoTracks.add("htracks_11_5_7_SelectedTrack_Axis_02_21", "htrack_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_5_8_SelectedTrack_Axis_03_21", "htrack_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_5_9_SelectedTrack_Axis_20_21", "htrack_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    recoTracks.add("htracks_11_6_1_SelectedTrack_Axis_00_22", "htrack_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // Axis_p             ;
    recoTracks.add("htracks_11_6_2_SelectedTrack_Axis_01_22", "htrack_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // Axis_pt            ;
    recoTracks.add("htracks_11_6_3_SelectedTrack_Axis_02_22", "htrack_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_6_4_SelectedTrack_Axis_03_22", "htrack_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_6_5_SelectedTrack_Axis_00_23", "htrack_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // Axis_p             ;
    recoTracks.add("htracks_11_6_6_SelectedTrack_Axis_01_23", "htrack_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // Axis_pt            ;
    recoTracks.add("htracks_11_6_7_SelectedTrack_Axis_02_23", "htrack_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_6_8_SelectedTrack_Axis_03_23", "htrack_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_6_9_SelectedTrack_Axis_22_23", "htrack_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    recoTracks.add("htracks_11_7_1_SelectedTrack_Axis_00_24", "htrack_Axis_00_24", kTH2F, {Axis_00, Axis_24}); // Axis_p             ;
    recoTracks.add("htracks_11_7_2_SelectedTrack_Axis_01_24", "htrack_Axis_01_24", kTH2F, {Axis_01, Axis_24}); // Axis_pt            ;
    recoTracks.add("htracks_11_7_3_SelectedTrack_Axis_02_24", "htrack_Axis_02_24", kTH2F, {Axis_02, Axis_24}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_7_4_SelectedTrack_Axis_03_24", "htrack_Axis_03_24", kTH2F, {Axis_03, Axis_24}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_7_5_SelectedTrack_Axis_00_25", "htrack_Axis_00_25", kTH2F, {Axis_00, Axis_25}); // Axis_p             ;
    recoTracks.add("htracks_11_7_6_SelectedTrack_Axis_01_25", "htrack_Axis_01_25", kTH2F, {Axis_01, Axis_25}); // Axis_pt            ;
    recoTracks.add("htracks_11_7_7_SelectedTrack_Axis_02_25", "htrack_Axis_02_25", kTH2F, {Axis_02, Axis_25}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_7_8_SelectedTrack_Axis_03_25", "htrack_Axis_03_25", kTH2F, {Axis_03, Axis_25}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_7_9_SelectedTrack_Axis_24_25", "htrack_Axis_24_25", kTH2F, {Axis_24, Axis_25}); // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    recoTracks.add("htracks_11_8_1_SelectedTrack_Axis_00_26", "htrack_Axis_00_26", kTH2F, {Axis_00, Axis_26}); // Axis_p             ;
    recoTracks.add("htracks_11_8_2_SelectedTrack_Axis_01_26", "htrack_Axis_01_26", kTH2F, {Axis_01, Axis_26}); // Axis_pt            ;
    recoTracks.add("htracks_11_8_3_SelectedTrack_Axis_02_26", "htrack_Axis_02_26", kTH2F, {Axis_02, Axis_26}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_8_4_SelectedTrack_Axis_03_26", "htrack_Axis_03_26", kTH2F, {Axis_03, Axis_26}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_8_5_SelectedTrack_Axis_00_27", "htrack_Axis_00_27", kTH2F, {Axis_00, Axis_27}); // Axis_p             ;
    recoTracks.add("htracks_11_8_6_SelectedTrack_Axis_01_27", "htrack_Axis_01_27", kTH2F, {Axis_01, Axis_27}); // Axis_pt            ;
    recoTracks.add("htracks_11_8_7_SelectedTrack_Axis_02_27", "htrack_Axis_02_27", kTH2F, {Axis_02, Axis_27}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_8_8_SelectedTrack_Axis_03_27", "htrack_Axis_03_27", kTH2F, {Axis_03, Axis_27}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_8_9_SelectedTrack_Axis_26_27", "htrack_Axis_26_27", kTH2F, {Axis_26, Axis_27}); // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    recoTracks.add("htracks_11_9_1_SelectedTrack_Axis_00_28", "htrack_Axis_00_28", kTH2F, {Axis_00, Axis_28}); // Axis_p             ;
    recoTracks.add("htracks_11_9_2_SelectedTrack_Axis_01_28", "htrack_Axis_01_28", kTH2F, {Axis_01, Axis_28}); // Axis_pt            ;
    recoTracks.add("htracks_11_9_3_SelectedTrack_Axis_02_28", "htrack_Axis_02_28", kTH2F, {Axis_02, Axis_28}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_9_4_SelectedTrack_Axis_03_28", "htrack_Axis_03_28", kTH2F, {Axis_03, Axis_28}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_9_5_SelectedTrack_Axis_00_29", "htrack_Axis_00_29", kTH2F, {Axis_00, Axis_29}); // Axis_p             ;
    recoTracks.add("htracks_11_9_6_SelectedTrack_Axis_01_29", "htrack_Axis_01_29", kTH2F, {Axis_01, Axis_29}); // Axis_pt            ;
    recoTracks.add("htracks_11_9_7_SelectedTrack_Axis_02_29", "htrack_Axis_02_29", kTH2F, {Axis_02, Axis_29}); // Axis_tpcInnerParam ;
    recoTracks.add("htracks_11_9_8_SelectedTrack_Axis_03_29", "htrack_Axis_03_29", kTH2F, {Axis_03, Axis_29}); // Axis_tofExpMom     ;
    recoTracks.add("htracks_11_9_9_SelectedTrack_Axis_28_29", "htrack_Axis_28_29", kTH2F, {Axis_28, Axis_29}); // Axis_tpcInnerParam ;
    // Deuteron
    // SelectedTrack
    //
    // Analysis
    // Pion
    // tpcSignal
    recoAnalysis.add("hAnalysis_11_Pi_Id0_1_Axis_00_05", "hAnalysis_11_Pi_Id0_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_1_Axis_02_05", "hAnalysis_11_Pi_Id0_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_1_Axis_03_05", "hAnalysis_11_Pi_Id0_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_11_Pi_Id0_2_Axis_00_06", "hAnalysis_11_Pi_Id0_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_2_Axis_02_06", "hAnalysis_11_Pi_Id0_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_2_Axis_03_06", "hAnalysis_11_Pi_Id0_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Pion
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_1_Axis_00_20", "hAnalysis_11_Pi_Id0_3_1_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // track.p            (),track.tpcNSigmaPi()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_2_Axis_01_20", "hAnalysis_11_Pi_Id0_3_2_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // track.pt           (),track.tpcNSigmaPi()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_3_Axis_02_20", "hAnalysis_11_Pi_Id0_3_3_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // track.tpcInnerParam(),track.tpcNSigmaPi()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_4_Axis_03_20", "hAnalysis_11_Pi_Id0_3_4_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // track.tofExpMom    (),track.tpcNSigmaPi()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_5_Axis_00_21", "hAnalysis_11_Pi_Id0_3_5_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // track.p            (),track.tofNSigmaPi()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_6_Axis_01_21", "hAnalysis_11_Pi_Id0_3_6_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // track.pt           (),track.tofNSigmaPi()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_7_Axis_02_21", "hAnalysis_11_Pi_Id0_3_7_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // track.tpcInnerParam(),track.tofNSigmaPi()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_8_Axis_03_21", "hAnalysis_11_Pi_Id0_3_8_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // track.tofExpMom    (),track.tofNSigmaPi()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_11_Pi_Id0_3_9_Axis_20_21", "hAnalysis_11_Pi_Id0_3_9_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // track.tpcNSigmaPi  (),track.tofNSigmaPi()) ;//Axis_tpcInnerParam ;

    // tpcSignal
    recoAnalysis.add("hAnalysis_11_Pi_Id1_1_Axis_00_05", "hAnalysis_11_Pi_Id1_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_1_Axis_02_05", "hAnalysis_11_Pi_Id1_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_1_Axis_03_05", "hAnalysis_11_Pi_Id1_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_11_Pi_Id1_2_Axis_00_06", "hAnalysis_11_Pi_Id1_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_2_Axis_02_06", "hAnalysis_11_Pi_Id1_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_2_Axis_03_06", "hAnalysis_11_Pi_Id1_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Pion
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_1_Axis_00_20", "hAnalysis_11_Pi_Id1_3_1_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // track.p            (),track.tpcNSigmaPi()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_2_Axis_01_20", "hAnalysis_11_Pi_Id1_3_2_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // track.pt           (),track.tpcNSigmaPi()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_3_Axis_02_20", "hAnalysis_11_Pi_Id1_3_3_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // track.tpcInnerParam(),track.tpcNSigmaPi()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_4_Axis_03_20", "hAnalysis_11_Pi_Id1_3_4_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // track.tofExpMom    (),track.tpcNSigmaPi()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_5_Axis_00_21", "hAnalysis_11_Pi_Id1_3_5_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // track.p            (),track.tofNSigmaPi()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_6_Axis_01_21", "hAnalysis_11_Pi_Id1_3_6_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // track.pt           (),track.tofNSigmaPi()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_7_Axis_02_21", "hAnalysis_11_Pi_Id1_3_7_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // track.tpcInnerParam(),track.tofNSigmaPi()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_8_Axis_03_21", "hAnalysis_11_Pi_Id1_3_8_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // track.tofExpMom    (),track.tofNSigmaPi()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_11_Pi_Id1_3_9_Axis_20_21", "hAnalysis_11_Pi_Id1_3_9_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // track.tpcNSigmaPi  (),track.tofNSigmaPi()) ;//Axis_tpcInnerParam ;
    // Pion

    // Kaon
    // tpcSignal
    recoAnalysis.add("hAnalysis_12_Ka_Id0_1_Axis_00_05", "hAnalysis_12_Ka_Id0_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_1_Axis_02_05", "hAnalysis_12_Ka_Id0_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_1_Axis_03_05", "hAnalysis_12_Ka_Id0_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_12_Ka_Id0_2_Axis_00_06", "hAnalysis_12_Ka_Id0_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_2_Axis_02_06", "hAnalysis_12_Ka_Id0_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_2_Axis_03_06", "hAnalysis_12_Ka_Id0_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Kaon
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_1_Axis_00_22", "hAnalysis_12_Ka_Id0_3_1_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // track.p            (),track.tpcNSigmaKa()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_2_Axis_01_22", "hAnalysis_12_Ka_Id0_3_2_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // track.pt           (),track.tpcNSigmaKa()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_3_Axis_02_22", "hAnalysis_12_Ka_Id0_3_3_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // track.tpcInnerParam(),track.tpcNSigmaKa()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_4_Axis_03_22", "hAnalysis_12_Ka_Id0_3_4_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // track.tofExpMom    (),track.tpcNSigmaKa()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_5_Axis_00_23", "hAnalysis_12_Ka_Id0_3_5_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // track.p            (),track.tofNSigmaKa()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_6_Axis_01_23", "hAnalysis_12_Ka_Id0_3_6_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // track.pt           (),track.tofNSigmaKa()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_7_Axis_02_23", "hAnalysis_12_Ka_Id0_3_7_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // track.tpcInnerParam(),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_8_Axis_03_23", "hAnalysis_12_Ka_Id0_3_8_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // track.tofExpMom    (),track.tofNSigmaKa()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_12_Ka_Id0_3_9_Axis_22_23", "hAnalysis_12_Ka_Id0_3_9_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // track.tpcNSigmaKa  (),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;

    // tpcSignal
    recoAnalysis.add("hAnalysis_12_Ka_Id1_1_Axis_00_05", "hAnalysis_12_Ka_Id1_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_1_Axis_02_05", "hAnalysis_12_Ka_Id1_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_1_Axis_03_05", "hAnalysis_12_Ka_Id1_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_12_Ka_Id1_2_Axis_00_06", "hAnalysis_12_Ka_Id1_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_2_Axis_02_06", "hAnalysis_12_Ka_Id1_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_2_Axis_03_06", "hAnalysis_12_Ka_Id1_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Kaon
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_1_Axis_00_22", "hAnalysis_12_Ka_Id1_3_1_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // track.p            (),track.tpcNSigmaKa()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_2_Axis_01_22", "hAnalysis_12_Ka_Id1_3_2_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // track.pt           (),track.tpcNSigmaKa()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_3_Axis_02_22", "hAnalysis_12_Ka_Id1_3_3_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // track.tpcInnerParam(),track.tpcNSigmaKa()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_4_Axis_03_22", "hAnalysis_12_Ka_Id1_3_4_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // track.tofExpMom    (),track.tpcNSigmaKa()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_5_Axis_00_23", "hAnalysis_12_Ka_Id1_3_5_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // track.p            (),track.tofNSigmaKa()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_6_Axis_01_23", "hAnalysis_12_Ka_Id1_3_6_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // track.pt           (),track.tofNSigmaKa()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_7_Axis_02_23", "hAnalysis_12_Ka_Id1_3_7_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // track.tpcInnerParam(),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_8_Axis_03_23", "hAnalysis_12_Ka_Id1_3_8_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // track.tofExpMom    (),track.tofNSigmaKa()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_12_Ka_Id1_3_9_Axis_22_23", "hAnalysis_12_Ka_Id1_3_9_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // track.tpcNSigmaKa  (),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;
    // Kaon

    // Proton
    // tpcSignal
    recoAnalysis.add("hAnalysis_13_Pr_Id0_1_Axis_00_05", "hAnalysis_13_Pr_Id0_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_1_Axis_02_05", "hAnalysis_13_Pr_Id0_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_1_Axis_03_05", "hAnalysis_13_Pr_Id0_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_13_Pr_Id0_2_Axis_00_06", "hAnalysis_13_Pr_Id0_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_2_Axis_02_06", "hAnalysis_13_Pr_Id0_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_2_Axis_03_06", "hAnalysis_13_Pr_Id0_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Proton
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_1_Axis_00_24", "hAnalysis_13_Pr_Id0_3_1_Axis_00_24", kTH2F, {Axis_00, Axis_24}); // track.p            (),track.tpcNSigmaPr()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_2_Axis_01_24", "hAnalysis_13_Pr_Id0_3_2_Axis_01_24", kTH2F, {Axis_01, Axis_24}); // track.pt           (),track.tpcNSigmaPr()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_3_Axis_02_24", "hAnalysis_13_Pr_Id0_3_3_Axis_02_24", kTH2F, {Axis_02, Axis_24}); // track.tpcInnerParam(),track.tpcNSigmaPr()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_4_Axis_03_24", "hAnalysis_13_Pr_Id0_3_4_Axis_03_24", kTH2F, {Axis_03, Axis_24}); // track.tofExpMom    (),track.tpcNSigmaPr()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_5_Axis_00_25", "hAnalysis_13_Pr_Id0_3_5_Axis_00_25", kTH2F, {Axis_00, Axis_25}); // track.p            (),track.tofNSigmaPr()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_6_Axis_01_25", "hAnalysis_13_Pr_Id0_3_6_Axis_01_25", kTH2F, {Axis_01, Axis_25}); // track.pt           (),track.tofNSigmaPr()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_7_Axis_02_25", "hAnalysis_13_Pr_Id0_3_7_Axis_02_25", kTH2F, {Axis_02, Axis_25}); // track.tpcInnerParam(),track.tofNSigmaPr()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_8_Axis_03_25", "hAnalysis_13_Pr_Id0_3_8_Axis_03_25", kTH2F, {Axis_03, Axis_25}); // track.tofExpMom    (),track.tofNSigmaPr()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_13_Pr_Id0_3_9_Axis_24_25", "hAnalysis_13_Pr_Id0_3_9_Axis_24_25", kTH2F, {Axis_24, Axis_25}); // track.tpcNSigmaPr  (),track.tofNSigmaPr()) ;//Axis_tpcInnerParam ;

    // tpcSignal
    recoAnalysis.add("hAnalysis_13_Pr_Id1_1_Axis_00_05", "hAnalysis_13_Pr_Id1_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_1_Axis_02_05", "hAnalysis_13_Pr_Id1_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_1_Axis_03_05", "hAnalysis_13_Pr_Id1_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_13_Pr_Id1_2_Axis_00_06", "hAnalysis_13_Pr_Id1_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_2_Axis_02_06", "hAnalysis_13_Pr_Id1_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_2_Axis_03_06", "hAnalysis_13_Pr_Id1_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Proton
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_1_Axis_00_24", "hAnalysis_13_Pr_Id1_3_1_Axis_00_24", kTH2F, {Axis_00, Axis_24}); // track.p            (),track.tpcNSigmaPr()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_2_Axis_01_24", "hAnalysis_13_Pr_Id1_3_2_Axis_01_24", kTH2F, {Axis_01, Axis_24}); // track.pt           (),track.tpcNSigmaPr()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_3_Axis_02_24", "hAnalysis_13_Pr_Id1_3_3_Axis_02_24", kTH2F, {Axis_02, Axis_24}); // track.tpcInnerParam(),track.tpcNSigmaPr()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_4_Axis_03_24", "hAnalysis_13_Pr_Id1_3_4_Axis_03_24", kTH2F, {Axis_03, Axis_24}); // track.tofExpMom    (),track.tpcNSigmaPr()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_5_Axis_00_25", "hAnalysis_13_Pr_Id1_3_5_Axis_00_25", kTH2F, {Axis_00, Axis_25}); // track.p            (),track.tofNSigmaPr()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_6_Axis_01_25", "hAnalysis_13_Pr_Id1_3_6_Axis_01_25", kTH2F, {Axis_01, Axis_25}); // track.pt           (),track.tofNSigmaPr()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_7_Axis_02_25", "hAnalysis_13_Pr_Id1_3_7_Axis_02_25", kTH2F, {Axis_02, Axis_25}); // track.tpcInnerParam(),track.tofNSigmaPr()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_8_Axis_03_25", "hAnalysis_13_Pr_Id1_3_8_Axis_03_25", kTH2F, {Axis_03, Axis_25}); // track.tofExpMom    (),track.tofNSigmaPr()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_13_Pr_Id1_3_9_Axis_24_25", "hAnalysis_13_Pr_Id1_3_9_Axis_24_25", kTH2F, {Axis_24, Axis_25}); // track.tpcNSigmaPr  (),track.tofNSigmaPr()) ;//Axis_tpcInnerParam ;
    // Proton

    // Electron
    // tpcSignal
    recoAnalysis.add("hAnalysis_14_El_Id0_1_Axis_00_05", "hAnalysis_14_El_Id0_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id0_1_Axis_02_05", "hAnalysis_14_El_Id0_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id0_1_Axis_03_05", "hAnalysis_14_El_Id0_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_14_El_Id0_2_Axis_00_06", "hAnalysis_14_El_Id0_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id0_2_Axis_02_06", "hAnalysis_14_El_Id0_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id0_2_Axis_03_06", "hAnalysis_14_El_Id0_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Electron
    recoAnalysis.add("hAnalysis_14_El_Id0_3_1_Axis_00_26", "hAnalysis_14_El_Id0_3_1_Axis_00_26", kTH2F, {Axis_00, Axis_26}); // track.p            (),track.tpcNSigmaEl()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_2_Axis_01_26", "hAnalysis_14_El_Id0_3_2_Axis_01_26", kTH2F, {Axis_01, Axis_26}); // track.pt           (),track.tpcNSigmaEl()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_3_Axis_02_26", "hAnalysis_14_El_Id0_3_3_Axis_02_26", kTH2F, {Axis_02, Axis_26}); // track.tpcInnerParam(),track.tpcNSigmaEl()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_4_Axis_03_26", "hAnalysis_14_El_Id0_3_4_Axis_03_26", kTH2F, {Axis_03, Axis_26}); // track.tofExpMom    (),track.tpcNSigmaEl()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_5_Axis_00_27", "hAnalysis_14_El_Id0_3_5_Axis_00_27", kTH2F, {Axis_00, Axis_27}); // track.p            (),track.tofNSigmaEl()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_6_Axis_01_27", "hAnalysis_14_El_Id0_3_6_Axis_01_27", kTH2F, {Axis_01, Axis_27}); // track.pt           (),track.tofNSigmaEl()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_7_Axis_02_27", "hAnalysis_14_El_Id0_3_7_Axis_02_27", kTH2F, {Axis_02, Axis_27}); // track.tpcInnerParam(),track.tofNSigmaEl()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_8_Axis_03_27", "hAnalysis_14_El_Id0_3_8_Axis_03_27", kTH2F, {Axis_03, Axis_27}); // track.tofExpMom    (),track.tofNSigmaEl()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_14_El_Id0_3_9_Axis_26_27", "hAnalysis_14_El_Id0_3_9_Axis_26_27", kTH2F, {Axis_26, Axis_27}); // track.tpcNSigmaEl  (),track.tofNSigmaEl()) ;//Axis_tpcInnerParam ;

    // tpcSignal
    recoAnalysis.add("hAnalysis_14_El_Id1_1_Axis_00_05", "hAnalysis_14_El_Id1_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id1_1_Axis_02_05", "hAnalysis_14_El_Id1_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id1_1_Axis_03_05", "hAnalysis_14_El_Id1_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_14_El_Id1_2_Axis_00_06", "hAnalysis_14_El_Id1_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id1_2_Axis_02_06", "hAnalysis_14_El_Id1_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id1_2_Axis_03_06", "hAnalysis_14_El_Id1_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Electron
    recoAnalysis.add("hAnalysis_14_El_Id1_3_1_Axis_00_26", "hAnalysis_14_El_Id1_3_1_Axis_00_26", kTH2F, {Axis_00, Axis_26}); // track.p            (),track.tpcNSigmaEl()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_2_Axis_01_26", "hAnalysis_14_El_Id1_3_2_Axis_01_26", kTH2F, {Axis_01, Axis_26}); // track.pt           (),track.tpcNSigmaEl()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_3_Axis_02_26", "hAnalysis_14_El_Id1_3_3_Axis_02_26", kTH2F, {Axis_02, Axis_26}); // track.tpcInnerParam(),track.tpcNSigmaEl()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_4_Axis_03_26", "hAnalysis_14_El_Id1_3_4_Axis_03_26", kTH2F, {Axis_03, Axis_26}); // track.tofExpMom    (),track.tpcNSigmaEl()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_5_Axis_00_27", "hAnalysis_14_El_Id1_3_5_Axis_00_27", kTH2F, {Axis_00, Axis_27}); // track.p            (),track.tofNSigmaEl()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_6_Axis_01_27", "hAnalysis_14_El_Id1_3_6_Axis_01_27", kTH2F, {Axis_01, Axis_27}); // track.pt           (),track.tofNSigmaEl()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_7_Axis_02_27", "hAnalysis_14_El_Id1_3_7_Axis_02_27", kTH2F, {Axis_02, Axis_27}); // track.tpcInnerParam(),track.tofNSigmaEl()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_8_Axis_03_27", "hAnalysis_14_El_Id1_3_8_Axis_03_27", kTH2F, {Axis_03, Axis_27}); // track.tofExpMom    (),track.tofNSigmaEl()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_14_El_Id1_3_9_Axis_26_27", "hAnalysis_14_El_Id1_3_9_Axis_26_27", kTH2F, {Axis_26, Axis_27}); // track.tpcNSigmaEl  (),track.tofNSigmaEl()) ;//Axis_tpcInnerParam ;
    // Electron

    // Deuteron
    // tpcSignal
    recoAnalysis.add("hAnalysis_15_De_Id0_1_Axis_00_05", "hAnalysis_15_De_Id0_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id0_1_Axis_02_05", "hAnalysis_15_De_Id0_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id0_1_Axis_03_05", "hAnalysis_15_De_Id0_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_15_De_Id0_2_Axis_00_06", "hAnalysis_15_De_Id0_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id0_2_Axis_02_06", "hAnalysis_15_De_Id0_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id0_2_Axis_03_06", "hAnalysis_15_De_Id0_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Deuteron
    recoAnalysis.add("hAnalysis_15_De_Id0_3_1_Axis_00_28", "hAnalysis_15_De_Id0_3_1_Axis_00_28", kTH2F, {Axis_00, Axis_28}); // track.p            (),track.tpcNSigmaDe()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_2_Axis_01_28", "hAnalysis_15_De_Id0_3_2_Axis_01_28", kTH2F, {Axis_01, Axis_28}); // track.pt           (),track.tpcNSigmaDe()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_3_Axis_02_28", "hAnalysis_15_De_Id0_3_3_Axis_02_28", kTH2F, {Axis_02, Axis_28}); // track.tpcInnerParam(),track.tpcNSigmaDe()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_4_Axis_03_28", "hAnalysis_15_De_Id0_3_4_Axis_03_28", kTH2F, {Axis_03, Axis_28}); // track.tofExpMom    (),track.tpcNSigmaDe()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_5_Axis_00_29", "hAnalysis_15_De_Id0_3_5_Axis_00_29", kTH2F, {Axis_00, Axis_29}); // track.p            (),track.tofNSigmaDe()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_6_Axis_01_29", "hAnalysis_15_De_Id0_3_6_Axis_01_29", kTH2F, {Axis_01, Axis_29}); // track.pt           (),track.tofNSigmaDe()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_7_Axis_02_29", "hAnalysis_15_De_Id0_3_7_Axis_02_29", kTH2F, {Axis_02, Axis_29}); // track.tpcInnerParam(),track.tofNSigmaDe()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_8_Axis_03_29", "hAnalysis_15_De_Id0_3_8_Axis_03_29", kTH2F, {Axis_03, Axis_29}); // track.tofExpMom    (),track.tofNSigmaDe()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_15_De_Id0_3_9_Axis_28_29", "hAnalysis_15_De_Id0_3_9_Axis_28_29", kTH2F, {Axis_28, Axis_29}); // track.tpcNSigmaDe  (),track.tofNSigmaDe()) ;//Axis_tpcInnerParam ;

    // tpcSignal
    recoAnalysis.add("hAnalysis_15_De_Id1_1_Axis_00_05", "hAnalysis_15_De_Id1_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id1_1_Axis_02_05", "hAnalysis_15_De_Id1_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id1_1_Axis_03_05", "hAnalysis_15_De_Id1_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    recoAnalysis.add("hAnalysis_15_De_Id1_2_Axis_00_06", "hAnalysis_15_De_Id1_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id1_2_Axis_02_06", "hAnalysis_15_De_Id1_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id1_2_Axis_03_06", "hAnalysis_15_De_Id1_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Deuteron
    recoAnalysis.add("hAnalysis_15_De_Id1_3_1_Axis_00_28", "hAnalysis_15_De_Id1_3_1_Axis_00_28", kTH2F, {Axis_00, Axis_28}); // track.p            (),track.tpcNSigmaDe()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_2_Axis_01_28", "hAnalysis_15_De_Id1_3_2_Axis_01_28", kTH2F, {Axis_01, Axis_28}); // track.pt           (),track.tpcNSigmaDe()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_3_Axis_02_28", "hAnalysis_15_De_Id1_3_3_Axis_02_28", kTH2F, {Axis_02, Axis_28}); // track.tpcInnerParam(),track.tpcNSigmaDe()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_4_Axis_03_28", "hAnalysis_15_De_Id1_3_4_Axis_03_28", kTH2F, {Axis_03, Axis_28}); // track.tofExpMom    (),track.tpcNSigmaDe()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_5_Axis_00_29", "hAnalysis_15_De_Id1_3_5_Axis_00_29", kTH2F, {Axis_00, Axis_29}); // track.p            (),track.tofNSigmaDe()) ;//Axis_p             ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_6_Axis_01_29", "hAnalysis_15_De_Id1_3_6_Axis_01_29", kTH2F, {Axis_01, Axis_29}); // track.pt           (),track.tofNSigmaDe()) ;//Axis_pt            ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_7_Axis_02_29", "hAnalysis_15_De_Id1_3_7_Axis_02_29", kTH2F, {Axis_02, Axis_29}); // track.tpcInnerParam(),track.tofNSigmaDe()) ;//Axis_tpcInnerParam ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_8_Axis_03_29", "hAnalysis_15_De_Id1_3_8_Axis_03_29", kTH2F, {Axis_03, Axis_29}); // track.tofExpMom    (),track.tofNSigmaDe()) ;//Axis_tofExpMom     ;
    recoAnalysis.add("hAnalysis_15_De_Id1_3_9_Axis_28_29", "hAnalysis_15_De_Id1_3_9_Axis_28_29", kTH2F, {Axis_28, Axis_29}); // track.tpcNSigmaDe  (),track.tofNSigmaDe()) ;//Axis_tpcInnerParam ;
    // Deuteron

    recoAnalysis.add("hAnalysis_4_SelectedTrack_IdentificationTag", "hAnalysis_SelectedTrack_IdentificationTag", {HistType::kTH1D, {{34, -1.5, 32.5, "trackTAG"}}});
    recoAnalysis.add("hAnalysis_4_RejectedTrack_RejectionTag", "hAnalysis_RejectedTrack_RejectionTag", {HistType::kTH1D, {{16, -1.5, 6.5, "rejectionTAG"}}});

    recoAnalysis.add("hAnalysis_5_Sparse_Full_K0sPiKa", "hAnalysis_5_Sparse_Full_K0sPiKa", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}});
    recoAnalysis.add("hAnalysis_5_Sparse_Full_K0sPrDe", "hAnalysis_5_Sparse_Full_K0sPrDe", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("hAnalysis_5_Sparse_Full_K0sKaEl", "hAnalysis_5_Sparse_Full_K0sKaEl", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}});
    recoAnalysis.add("hAnalysis_5_Sparse_Full_PiKaPr", "hAnalysis_5_Sparse_Full_PiKaPr", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}});
    recoAnalysis.add("hAnalysis_5_Sparse_Full_PiElDe", "hAnalysis_5_Sparse_Full_PiElDe", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("hAnalysis_5_Sparse_Full_KaPrDe", "hAnalysis_5_Sparse_Full_KaPrDe", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("hAnalysis_5_Sparse_Full_PrElDe", "hAnalysis_5_Sparse_Full_PrElDe", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});

    recoAnalysis.add("hAnalysis_5_Sparse_newDynm_K0s_Ka", "hAnalysis_Sparse_newDynm_K0s_Ka", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {500, -1.5, 498.5, "nKaon"}, {10000, -1.5, 9998.5, "(nK0s)^{2}"}, {250000, -1.5, 249998.5, "(nKaon)^{2}"}, {500, -1.5, 498.5, "(nK0s*nKaon)"}});
    recoAnalysis.add("hAnalysis_5_Sparse_newDynm_Kp_Km", "hAnalysis_Sparse_newDynm_Kp_Km", kTHnSparseD, {Axis_centFT0C, {2000, -1.5, 1998.5, "nTrack"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {250000, -1.5, 249998.5, "(nKaPlus)^{2}"}, {250000, -1.5, 249998.5, "(nKaMinus)^{2}"}, {250000, -1.5, 249998.5, "(nKaPlus*nKaMinus)"}});
    recoAnalysis.add("hTest", "hTest", kTH1D, {VectorAxis});
    //
    // Printing the Stored Registry information
    LOG(info) << "Printing Stored Registry Information";
    LOG(info) << "Printing RecoEvent ";
    recoEvent.print();
    LOG(info) << "Printing RecoK0s ";
    recoK0s.print();
    LOG(info) << "Printing RecoAnalysis ";
    recoAnalysis.print();
  }

  // tpc Selections
  template <typename T>
  bool selPionTPCInnerParam(T track)
  {
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaKa()) > 3.0 && TMath::Abs(track.tpcNSigmaPr()) > 3.0 && TMath::Abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && TMath::Abs(track.tpcNSigmaPi()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaPi()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selKaonTPCInnerParam(T track)
  {
    // p dependent cuts
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaPi()) > 3.0 && TMath::Abs(track.tpcNSigmaPr()) > 3.0 && TMath::Abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && TMath::Abs(track.tpcNSigmaKa()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaKa()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selProtonTPCInnerParam(T track)
  {
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaPi()) > 3.0 && TMath::Abs(track.tpcNSigmaKa()) > 3.0 && TMath::Abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.60 && TMath::Abs(track.tpcNSigmaPr()) < 3.0) {
        return true;
      }
      if (1.60 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaPr()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selDeuteronTPCInnerParam(T track)
  {
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaPi()) > 3.0 && TMath::Abs(track.tpcNSigmaKa()) > 3.0 && TMath::Abs(track.tpcNSigmaPr()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.80 && TMath::Abs(track.tpcNSigmaDe()) < 3.0) {
        return true;
      }
      if (1.80 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaDe()) < 2.0) {
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
    if (track.p() <= 0.75
        // && (TMath::Power(track.tpcNSigmaPi(),2)+TMath::Power(track.tofNSigmaPi(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaPi()) < 3.0 && TMath::Abs(track.tofNSigmaPi()) < 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    } else if (0.75 < track.p() // after p = 0.75, Pi and Ka lines of nSigma 3.0 will start intersecting
                                //  && (TMath::Power(track.tpcNSigmaPi(),2)+TMath::Power(track.tofNSigmaPi(),2)) < 8.0
               && TMath::Abs(track.tpcNSigmaPi()) < 2.0 && TMath::Abs(track.tofNSigmaPi()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaonTOF(T track)
  {
    if (track.p() <= 0.75
        // && (TMath::Power(track.tpcNSigmaKa(),2)+TMath::Power(track.tofNSigmaKa(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaKa()) < 3.0 && TMath::Abs(track.tofNSigmaKa()) < 3.0) {
      return true;
    }
    if (0.75 < track.p() && track.p() <= 1.30 // after 0.75 Pi and Ka lines of nSigma 3.0 will start intersecting
        && TMath::Abs(track.tpcNSigmaKa()) < 3.0 && TMath::Abs(track.tofNSigmaKa()) < 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    }
    if (1.30 < track.p() // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
                         //  && (TMath::Power(track.tpcNSigmaKa(),2)+TMath::Power(track.tofNSigmaKa(),2)) < 8.0
        && TMath::Abs(track.tpcNSigmaKa()) < 2.0 && TMath::Abs(track.tofNSigmaKa()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProtonTOF(T track)
  {
    if (track.p() <= 1.30
        // && (TMath::Power(track.tpcNSigmaPr(),2)+TMath::Power(track.tofNSigmaPr(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaPr()) < 3.0 && TMath::Abs(track.tofNSigmaPr()) < 3.0) {
      return true;
    }
    if (1.30 < track.p() && track.p() <= 3.10                                                                                                                                                                                                                 // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
        && TMath::Abs(track.tpcNSigmaPr()) < 3.0 && TMath::Abs(track.tofNSigmaPr()) < 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0 // Some Deuteron contamination is still coming in p dependent cuts
    ) {
      return true;
    }
    if (3.10 < track.p() // after 3.10 Pr and De lines of nSigma 3.0 will start intersecting
                         //  && (TMath::Power(track.tpcNSigmaPr(),2)+TMath::Power(track.tofNSigmaPr(),2)) < 8.0
        && TMath::Abs(track.tpcNSigmaPr()) < 2.0 && TMath::Abs(track.tofNSigmaPr()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteronTOF(T track)
  {
    if (track.p() <= 3.10
        // && (TMath::Power(track.tpcNSigmaDe(),2)+TMath::Power(track.tofNSigmaDe(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaDe()) < 3.0 && TMath::Abs(track.tofNSigmaDe()) < 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0) {
      return true;
    }
    if (3.10 < track.p()                                                                                                                                                                                                                                      // after 3.10 De and Pr lines of nSigma 3.0 will start intersecting
        && TMath::Abs(track.tpcNSigmaDe()) < 2.0 && TMath::Abs(track.tofNSigmaDe()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 // Some Deuteron contamination is still coming in p dependent cuts
    ) {
      return true;
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectronTOF(T track)
  {
    if (
      (TMath::Power(track.tpcNSigmaEl(), 2) + TMath::Power(track.tofNSigmaEl(), 2)) < 9.00 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
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
  //

  enum pidTagValue {
    tagOther = 0,
    tagPion = 1,
    tagKaon = 2,
    tagProton = 4,
    tagElectron = 8,
    tagDeuteron = 16
  };

  template <typename T>
  int FindTrackTag(T track)
  {
    int tempPID = tagOther;
    // 0- Some other particle
    // 1- Pion
    // 2- Kaon
    // 4- Proton
    // 8- Electron
    int PiIdMethod = -1;
    int KaIdMethod = -1;
    int PrIdMethod = -1;
    int ElIdMethod = -1;
    int DeIdMethod = -1;
    if (selPion(track, PiIdMethod)) {
      tempPID += tagPion;
    }
    if (selKaon(track, KaIdMethod)) {
      tempPID += tagKaon;
    }
    if (selProton(track, PrIdMethod)) {
      tempPID += tagProton;
    }
    if (selElectron(track, ElIdMethod)) {
      tempPID += tagElectron;
    }
    if (selDeuteron(track, DeIdMethod)) {
      tempPID += tagDeuteron;
    }

    return tempPID;
  }

  // V0 PID checks
  template <typename T>
  bool V0SelPion(T track, int& IdMethod)
  {
    if (selPion(track, IdMethod)) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool V0SelProton(T track, int& IdMethod)
  {
    if (selProton(track, IdMethod)) {
      return true;
    }
    return false;
  }

  template <typename T>
  int FindV0Tag(T posDaughterTrack, T negDaughterTrack, int& posPiIdMethod, int& posPrIdMethod, int& negPiIdMethod, int& negPrIdMethod)
  {

    bool posIsPion = false;
    bool posIsProton = false;
    bool negIsPion = false;
    bool negIsProton = false;

    int V0TagValue = 0;

    // Check if positive track is pion or proton
    if (V0SelPion(posDaughterTrack, posPiIdMethod))
      posIsPion = true; // Coming From K0s    -> PiPlus + PiMinus and AntiLambda -> PiPlus + AntiProton
    if (V0SelProton(posDaughterTrack, posPrIdMethod))
      posIsProton = true; // Coming From Lambda -> proton + PiMinus
    if (V0SelPion(negDaughterTrack, negPiIdMethod))
      negIsPion = true; // Coming From K0s       -> PiPlus + PiMinus and Lambda -> proton + PiMinus
    if (V0SelProton(negDaughterTrack, negPrIdMethod))
      negIsProton = true; // Coming From AntiLambda -> PiPlus + AntiProton

    if (posIsPion && negIsPion) {
      V0TagValue += 1;
    } // It is K0s
    if (posIsProton && negIsPion) {
      V0TagValue += 2;
    } // It is Lambda
    if (posIsPion && negIsProton) {
      V0TagValue += 4;
    } // It is AntiLambda

    return V0TagValue;
  }

  template <typename T>
  bool selK0s(T v0)
  {
    if (mLowK0s < v0.mK0Short() && v0.mK0Short() < mHighK0s && 0.1 < v0.pt() && v0.pt() < 1.5
        //  && TMath::Abs(v0.eta()) < 0.5
        && TMath::Abs(v0.rapidity(pdgDB->Mass(310))) < 0.5) {
      return true;
    } else {
      return false;
    }
  }

  // Insertion Sort  // Using insertion sort when array is almost sorted
  void InsertionSortVector(std::vector<int64_t>& UnsortedVector)
  {
    for (long unsigned int i = 1; i < UnsortedVector.size(); i++) {
      int currentElement = UnsortedVector[i]; // Element to be Inserted at correct position
      int j;                                  //(j+1) is the correct position of current element
      for (j = i - 1; j >= 0 && (UnsortedVector[j] > currentElement); j--) {
        UnsortedVector[j + 1] = UnsortedVector[j];
      }
      UnsortedVector[j + 1] = currentElement;
    }
  }

  // 03-Sorting a Vector Of Vector
  void InsertionSortVectorOfVector(std::vector<std::vector<int64_t>>& UnsortedVector)
  {
    for (long unsigned int i = 1; i < UnsortedVector.size(); i++) {
      std::vector<int64_t> currentElement = UnsortedVector[i]; // Element to be Inserted at correct position
      int j;                                                   //(j+1) is the correct position of current element
      for (j = i - 1; j >= 0 && (UnsortedVector[j][0] > currentElement[0]); j--) {
        UnsortedVector[j + 1] = UnsortedVector[j];
      }
      UnsortedVector[j + 1] = currentElement;
    }
  }

  template <typename T>
  long int FindIteratorPosition(long int globalIDX, T tracks)
  {
    long int iterPos = 0;
    for (auto trk : tracks) {
      if (trk.globalIndex() == globalIDX) {
        return iterPos;
      } else if (trk.globalIndex() > globalIDX) {
        return -1;
      }
      iterPos++;
    }
    return -1; // -1 means track is not in the table
  }

  template <typename T, typename U, typename V>
  void fillV0_K0sFullInformation(T v0, U posDaughterTrack, V negDaughterTrack, int posIdMethod, int negIdMethod)
  {
    recoV0s.fill(HIST("hV0s_table_0_2_01_K0s_Mass_Selected"), v0.mK0Short());
    recoV0s.fill(HIST("hV0s_table_0_2_02_K0s_P"), v0.p());
    recoV0s.fill(HIST("hV0s_table_0_2_03_K0s_Pt"), v0.pt());
    recoV0s.fill(HIST("hV0s_table_0_2_04_K0s_Eta"), v0.eta());
    recoV0s.fill(HIST("hV0s_table_0_2_05_K0s_Phi"), v0.phi());
    recoV0s.fill(HIST("hV0s_table_0_2_06_K0s_Rapidity"), v0.rapidity(pdgDB->Mass(310)));

    // posDaughterTrack
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_P"), posDaughterTrack.p());
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_Pt"), posDaughterTrack.pt());                       // Axis_pt            ;
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_tpcInnerParam"), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_tofExpMom"), posDaughterTrack.tofExpMom());         // Axis_tofExpMom     ;

    recoV0s.fill(HIST("hV0s_table_0_3_02_K0s_Pi_Eta"), posDaughterTrack.eta());
    recoV0s.fill(HIST("hV0s_table_0_3_03_K0s_Pi_Phi"), posDaughterTrack.phi());
    recoV0s.fill(HIST("hV0s_table_0_3_04_K0s_Pi_Rapidity"), posDaughterTrack.rapidity(pdgDB->Mass(211)));
    recoV0s.fill(HIST("hV0s_table_0_3_05_K0s_Pi_isPVContributor"), posDaughterTrack.isPVContributor());
    recoV0s.fill(HIST("hV0s_table_0_3_06_K0s_Pi_isGlobalTrack"), posDaughterTrack.isGlobalTrack());
    recoV0s.fill(HIST("hV0s_table_0_3_07_K0s_Pi_DcaXY"), posDaughterTrack.dcaXY());
    recoV0s.fill(HIST("hV0s_table_0_3_08_K0s_Pi_DcaZ"), posDaughterTrack.dcaZ());
    if (posIdMethod == 0) {
      // momemtum
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_01"), posDaughterTrack.p(), posDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_02"), posDaughterTrack.p(), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_03"), posDaughterTrack.p(), posDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_05"), posDaughterTrack.p(), posDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_05"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_05"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_06"), posDaughterTrack.p(), posDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_06"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_06"), posDaughterTrack.tofExpMom(), posDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_20"), posDaughterTrack.p(), posDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_20"), posDaughterTrack.pt(), posDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_20"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_20"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_21"), posDaughterTrack.p(), posDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_21"), posDaughterTrack.pt(), posDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_21"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_21"), posDaughterTrack.tofExpMom(), posDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_20_21"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    if (posIdMethod == 1) {
      // momemtum
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_01"), posDaughterTrack.p(), posDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_02"), posDaughterTrack.p(), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_03"), posDaughterTrack.p(), posDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_05"), posDaughterTrack.p(), posDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_05"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_05"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_06"), posDaughterTrack.p(), posDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_06"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_06"), posDaughterTrack.tofExpMom(), posDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_20"), posDaughterTrack.p(), posDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_20"), posDaughterTrack.pt(), posDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_20"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_20"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_21"), posDaughterTrack.p(), posDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_21"), posDaughterTrack.pt(), posDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_21"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_21"), posDaughterTrack.tofExpMom(), posDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_20_21"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    //
    // negDaughterTrack
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_P"), negDaughterTrack.p());
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_Pt"), negDaughterTrack.pt());                       // Axis_pt            ;
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_tpcInnerParam"), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoV0s.fill(HIST("hV0s_table_0_3_01_K0s_Pi_tofExpMom"), negDaughterTrack.tofExpMom());         // Axis_tofExpMom     ;

    recoV0s.fill(HIST("hV0s_table_0_3_02_K0s_Pi_Eta"), negDaughterTrack.eta());
    recoV0s.fill(HIST("hV0s_table_0_3_03_K0s_Pi_Phi"), negDaughterTrack.phi());
    recoV0s.fill(HIST("hV0s_table_0_3_04_K0s_Pi_Rapidity"), negDaughterTrack.rapidity(pdgDB->Mass(211)));
    recoV0s.fill(HIST("hV0s_table_0_3_05_K0s_Pi_isPVContributor"), negDaughterTrack.isPVContributor());
    recoV0s.fill(HIST("hV0s_table_0_3_06_K0s_Pi_isGlobalTrack"), negDaughterTrack.isGlobalTrack());
    recoV0s.fill(HIST("hV0s_table_0_3_07_K0s_Pi_DcaXY"), negDaughterTrack.dcaXY());
    recoV0s.fill(HIST("hV0s_table_0_3_08_K0s_Pi_DcaZ"), negDaughterTrack.dcaZ());
    if (negIdMethod == 0) {
      // momemtum
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_01"), negDaughterTrack.p(), negDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_02"), negDaughterTrack.p(), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_03"), negDaughterTrack.p(), negDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_05"), negDaughterTrack.p(), negDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_05"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_05"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_06"), negDaughterTrack.p(), negDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_06"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_06"), negDaughterTrack.tofExpMom(), negDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_20"), negDaughterTrack.p(), negDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_20"), negDaughterTrack.pt(), negDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_20"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_20"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_00_21"), negDaughterTrack.p(), negDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_01_21"), negDaughterTrack.pt(), negDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_02_21"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_03_21"), negDaughterTrack.tofExpMom(), negDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id0_Axis_20_21"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    if (negIdMethod == 1) {
      // momemtum
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_01"), negDaughterTrack.p(), negDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_02"), negDaughterTrack.p(), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_03"), negDaughterTrack.p(), negDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_05"), negDaughterTrack.p(), negDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_05"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_05"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_06"), negDaughterTrack.p(), negDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_06"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_06"), negDaughterTrack.tofExpMom(), negDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_20"), negDaughterTrack.p(), negDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_20"), negDaughterTrack.pt(), negDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_20"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_20"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_00_21"), negDaughterTrack.p(), negDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_01_21"), negDaughterTrack.pt(), negDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_02_21"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_03_21"), negDaughterTrack.tofExpMom(), negDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_table_0_3_10_K0s_Pi_Id1_Axis_20_21"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    //
  }

  template <typename T, typename U, typename V>
  void fillV0_K0sSelectedInformation(T v0, U posDaughterTrack, V negDaughterTrack, int posIdMethod, int negIdMethod)
  {
    recoV0s.fill(HIST("hV0s_0_2_01_K0s_Mass_Selected_LambdaMass"), v0.mLambda());
    recoV0s.fill(HIST("hV0s_0_2_01_K0s_Mass_Selected_AntiLambdaMass"), v0.mAntiLambda());
    recoV0s.fill(HIST("hV0s_0_2_01_K0s_Mass_Selected"), v0.mK0Short());
    recoV0s.fill(HIST("hV0s_0_2_02_K0s_P"), v0.p());
    recoV0s.fill(HIST("hV0s_0_2_03_K0s_Pt"), v0.pt());
    recoV0s.fill(HIST("hV0s_0_2_04_K0s_Eta"), v0.eta());
    recoV0s.fill(HIST("hV0s_0_2_05_K0s_Phi"), v0.phi());
    recoV0s.fill(HIST("hV0s_0_2_06_K0s_Rapidity"), v0.rapidity(pdgDB->Mass(310)));

    // posDaughterTrack
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_P"), posDaughterTrack.p());
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_Pt"), posDaughterTrack.pt());                       // Axis_pt            ;
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_tpcInnerParam"), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_tofExpMom"), posDaughterTrack.tofExpMom());         // Axis_tofExpMom     ;

    recoV0s.fill(HIST("hV0s_0_3_02_K0s_Pi_Eta"), posDaughterTrack.eta());
    recoV0s.fill(HIST("hV0s_0_3_03_K0s_Pi_Phi"), posDaughterTrack.phi());
    recoV0s.fill(HIST("hV0s_0_3_04_K0s_Pi_Rapidity"), posDaughterTrack.rapidity(pdgDB->Mass(211)));
    recoV0s.fill(HIST("hV0s_0_3_05_K0s_Pi_isPVContributor"), posDaughterTrack.isPVContributor());
    recoV0s.fill(HIST("hV0s_0_3_06_K0s_Pi_isGlobalTrack"), posDaughterTrack.isGlobalTrack());
    recoV0s.fill(HIST("hV0s_0_3_07_K0s_Pi_DcaXY"), posDaughterTrack.dcaXY());
    recoV0s.fill(HIST("hV0s_0_3_08_K0s_Pi_DcaZ"), posDaughterTrack.dcaZ());
    if (posIdMethod == 0) {
      // momemtum
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_01"), posDaughterTrack.p(), posDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_02"), posDaughterTrack.p(), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_03"), posDaughterTrack.p(), posDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_05"), posDaughterTrack.p(), posDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_05"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_05"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_06"), posDaughterTrack.p(), posDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_06"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_06"), posDaughterTrack.tofExpMom(), posDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_20"), posDaughterTrack.p(), posDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_01_20"), posDaughterTrack.pt(), posDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_20"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_20"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_21"), posDaughterTrack.p(), posDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_01_21"), posDaughterTrack.pt(), posDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_21"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_21"), posDaughterTrack.tofExpMom(), posDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_20_21"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    if (posIdMethod == 1) {
      // momemtum
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_01"), posDaughterTrack.p(), posDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_02"), posDaughterTrack.p(), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_03"), posDaughterTrack.p(), posDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_05"), posDaughterTrack.p(), posDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_05"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_05"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_06"), posDaughterTrack.p(), posDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_06"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_06"), posDaughterTrack.tofExpMom(), posDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_20"), posDaughterTrack.p(), posDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_01_20"), posDaughterTrack.pt(), posDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_20"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_20"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_21"), posDaughterTrack.p(), posDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_01_21"), posDaughterTrack.pt(), posDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_21"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_21"), posDaughterTrack.tofExpMom(), posDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_20_21"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    //

    // negDaughterTrack
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_P"), negDaughterTrack.p());
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_Pt"), negDaughterTrack.pt());                       // Axis_pt            ;
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_tpcInnerParam"), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoV0s.fill(HIST("hV0s_0_3_01_K0s_Pi_tofExpMom"), negDaughterTrack.tofExpMom());         // Axis_tofExpMom     ;

    recoV0s.fill(HIST("hV0s_0_3_02_K0s_Pi_Eta"), negDaughterTrack.eta());
    recoV0s.fill(HIST("hV0s_0_3_03_K0s_Pi_Phi"), negDaughterTrack.phi());
    recoV0s.fill(HIST("hV0s_0_3_04_K0s_Pi_Rapidity"), negDaughterTrack.rapidity(pdgDB->Mass(211)));
    recoV0s.fill(HIST("hV0s_0_3_05_K0s_Pi_isPVContributor"), negDaughterTrack.isPVContributor());
    recoV0s.fill(HIST("hV0s_0_3_06_K0s_Pi_isGlobalTrack"), negDaughterTrack.isGlobalTrack());
    recoV0s.fill(HIST("hV0s_0_3_07_K0s_Pi_DcaXY"), negDaughterTrack.dcaXY());
    recoV0s.fill(HIST("hV0s_0_3_08_K0s_Pi_DcaZ"), negDaughterTrack.dcaZ());
    if (negIdMethod == 0) {
      // momemtum
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_01"), negDaughterTrack.p(), negDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_02"), negDaughterTrack.p(), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_03"), negDaughterTrack.p(), negDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_05"), negDaughterTrack.p(), negDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_05"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_05"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_06"), negDaughterTrack.p(), negDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_06"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_06"), negDaughterTrack.tofExpMom(), negDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_20"), negDaughterTrack.p(), negDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_01_20"), negDaughterTrack.pt(), negDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_20"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_20"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_00_21"), negDaughterTrack.p(), negDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_01_21"), negDaughterTrack.pt(), negDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_02_21"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_03_21"), negDaughterTrack.tofExpMom(), negDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id0_Axis_20_21"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    if (negIdMethod == 1) {
      // momemtum
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_01"), negDaughterTrack.p(), negDaughterTrack.pt());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_02"), negDaughterTrack.p(), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_03"), negDaughterTrack.p(), negDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_05"), negDaughterTrack.p(), negDaughterTrack.tpcSignal());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_05"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_05"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_06"), negDaughterTrack.p(), negDaughterTrack.beta());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_06"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_06"), negDaughterTrack.tofExpMom(), negDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_20"), negDaughterTrack.p(), negDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_01_20"), negDaughterTrack.pt(), negDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_20"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_20"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_00_21"), negDaughterTrack.p(), negDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_01_21"), negDaughterTrack.pt(), negDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_02_21"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_03_21"), negDaughterTrack.tofExpMom(), negDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoV0s.fill(HIST("hV0s_0_3_10_K0s_Pi_Id1_Axis_20_21"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    //

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
    recoV0s.fill(HIST("hV0s_0_4_01_K0s_v0DaughterCollisionIndexTag"), v0daughterCollisionIndexTag);
  }

  template <typename T, typename U, typename V>
  void fillK0sSelectedInformation(T v0, U posDaughterTrack, V negDaughterTrack, int posIdMethod, int negIdMethod)
  {
    recoK0s.fill(HIST("hK0s_2_01_K0s_Mass_Selected_LambdaMass"), v0.mLambda());
    recoK0s.fill(HIST("hK0s_2_01_K0s_Mass_Selected_AntiLambdaMass"), v0.mAntiLambda());
    recoK0s.fill(HIST("hK0s_2_01_Mass_Selected"), v0.mK0Short());
    recoK0s.fill(HIST("hK0s_2_02_P"), v0.p());
    recoK0s.fill(HIST("hK0s_2_03_Pt"), v0.pt());
    recoK0s.fill(HIST("hK0s_2_04_Eta"), v0.eta());
    recoK0s.fill(HIST("hK0s_2_05_Phi"), v0.phi());
    recoK0s.fill(HIST("hK0s_2_06_Rapidity"), v0.rapidity(pdgDB->Mass(310)));

    // posDaughterTrack
    recoK0s.fill(HIST("hK0s_3_01_Pion_P"), posDaughterTrack.p());
    recoK0s.fill(HIST("hK0s_3_01_Pion_Pt"), posDaughterTrack.pt());                       // Axis_pt            ;
    recoK0s.fill(HIST("hK0s_3_01_Pion_tpcInnerParam"), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoK0s.fill(HIST("hK0s_3_01_Pion_tofExpMom"), posDaughterTrack.tofExpMom());         // Axis_tofExpMom     ;

    recoK0s.fill(HIST("hK0s_3_02_Pion_Eta"), posDaughterTrack.eta());
    recoK0s.fill(HIST("hK0s_3_03_Pion_Phi"), posDaughterTrack.phi());
    recoK0s.fill(HIST("hK0s_3_04_Pion_Rapidity"), posDaughterTrack.rapidity(pdgDB->Mass(211)));
    recoK0s.fill(HIST("hK0s_3_05_Daughter_isPVContributor"), posDaughterTrack.isPVContributor());
    recoK0s.fill(HIST("hK0s_3_06_Daughter_isGlobalTrack"), posDaughterTrack.isGlobalTrack());
    recoK0s.fill(HIST("hK0s_3_07_Daughter_DcaXY"), posDaughterTrack.dcaXY());
    recoK0s.fill(HIST("hK0s_3_08_Daughter_DcaZ"), posDaughterTrack.dcaZ());

    if (posIdMethod == 0) {
      // momemtum
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_01"), posDaughterTrack.p(), posDaughterTrack.pt());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_02"), posDaughterTrack.p(), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_03"), posDaughterTrack.p(), posDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_05"), posDaughterTrack.p(), posDaughterTrack.tpcSignal());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_05"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_05"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_06"), posDaughterTrack.p(), posDaughterTrack.beta());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_06"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_06"), posDaughterTrack.tofExpMom(), posDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_20"), posDaughterTrack.p(), posDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_01_20"), posDaughterTrack.pt(), posDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_20"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_20"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_21"), posDaughterTrack.p(), posDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_01_21"), posDaughterTrack.pt(), posDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_21"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_21"), posDaughterTrack.tofExpMom(), posDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_20_21"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    if (posIdMethod == 1) {
      // momemtum
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_01"), posDaughterTrack.p(), posDaughterTrack.pt());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_02"), posDaughterTrack.p(), posDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_03"), posDaughterTrack.p(), posDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_05"), posDaughterTrack.p(), posDaughterTrack.tpcSignal());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_05"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_05"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_06"), posDaughterTrack.p(), posDaughterTrack.beta());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_06"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_06"), posDaughterTrack.tofExpMom(), posDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_20"), posDaughterTrack.p(), posDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_01_20"), posDaughterTrack.pt(), posDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_20"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_20"), posDaughterTrack.tofExpMom(), posDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_21"), posDaughterTrack.p(), posDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_01_21"), posDaughterTrack.pt(), posDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_21"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_21"), posDaughterTrack.tofExpMom(), posDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_20_21"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    // negDaughterTrack
    recoK0s.fill(HIST("hK0s_3_01_Pion_P"), negDaughterTrack.p());
    recoK0s.fill(HIST("hK0s_3_01_Pion_Pt"), negDaughterTrack.pt());                       // Axis_pt            ;
    recoK0s.fill(HIST("hK0s_3_01_Pion_tpcInnerParam"), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoK0s.fill(HIST("hK0s_3_01_Pion_tofExpMom"), negDaughterTrack.tofExpMom());         // Axis_tofExpMom     ;

    recoK0s.fill(HIST("hK0s_3_02_Pion_Eta"), negDaughterTrack.eta());
    recoK0s.fill(HIST("hK0s_3_03_Pion_Phi"), negDaughterTrack.phi());
    recoK0s.fill(HIST("hK0s_3_04_Pion_Rapidity"), negDaughterTrack.rapidity(pdgDB->Mass(211)));
    recoK0s.fill(HIST("hK0s_3_05_Daughter_isPVContributor"), negDaughterTrack.isPVContributor());
    recoK0s.fill(HIST("hK0s_3_06_Daughter_isGlobalTrack"), negDaughterTrack.isGlobalTrack());
    recoK0s.fill(HIST("hK0s_3_07_Daughter_DcaXY"), negDaughterTrack.dcaXY());
    recoK0s.fill(HIST("hK0s_3_08_Daughter_DcaZ"), negDaughterTrack.dcaZ());

    if (negIdMethod == 0) {
      // momemtum
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_01"), negDaughterTrack.p(), negDaughterTrack.pt());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_02"), negDaughterTrack.p(), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_03"), negDaughterTrack.p(), negDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_05"), negDaughterTrack.p(), negDaughterTrack.tpcSignal());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_05"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_05"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_06"), negDaughterTrack.p(), negDaughterTrack.beta());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_06"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_06"), negDaughterTrack.tofExpMom(), negDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_20"), negDaughterTrack.p(), negDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_01_20"), negDaughterTrack.pt(), negDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_20"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_20"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_00_21"), negDaughterTrack.p(), negDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_01_21"), negDaughterTrack.pt(), negDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_02_21"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_03_21"), negDaughterTrack.tofExpMom(), negDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id0_Axis_20_21"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    if (negIdMethod == 1) {
      // momemtum
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_01"), negDaughterTrack.p(), negDaughterTrack.pt());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_02"), negDaughterTrack.p(), negDaughterTrack.tpcInnerParam()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_03"), negDaughterTrack.p(), negDaughterTrack.tofExpMom());     // Axis_tofExpMom     ;
      // tpcSignal
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_05"), negDaughterTrack.p(), negDaughterTrack.tpcSignal());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_05"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcSignal()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_05"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_06"), negDaughterTrack.p(), negDaughterTrack.beta());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_06"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.beta()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_06"), negDaughterTrack.tofExpMom(), negDaughterTrack.beta());     // Axis_tofExpMom     ;
      // Look at Pion
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_20"), negDaughterTrack.p(), negDaughterTrack.tpcNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_01_20"), negDaughterTrack.pt(), negDaughterTrack.tpcNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_20"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_20"), negDaughterTrack.tofExpMom(), negDaughterTrack.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_00_21"), negDaughterTrack.p(), negDaughterTrack.tofNSigmaPi());             // Axis_p             ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_01_21"), negDaughterTrack.pt(), negDaughterTrack.tofNSigmaPi());            // Axis_pt            ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_02_21"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_03_21"), negDaughterTrack.tofExpMom(), negDaughterTrack.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoK0s.fill(HIST("hK0s_3_10_K0s_Pi_Id1_Axis_20_21"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
    //

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

    recoK0s.fill(HIST("hK0s_4_01_v0DaughterCollisionIndexTag"), v0daughterCollisionIndexTag);
  }

  template <typename T>
  void FindRepeatEntries(std::vector<int64_t> ParticleList, T hist)
  {
    for (long unsigned int ii = 0; ii < ParticleList.size(); ii++) {
      int nCommonCount = 0; // checking the repeat number of track
      for (long unsigned int jj = 0; jj < ParticleList.size(); jj++) {
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

  void FillNewListFromOldList(std::vector<int64_t>& NewList, std::vector<int64_t> OldList)
  {
    for (long unsigned int ii = 0; ii < OldList.size(); ii++) {
      bool RepeatEntry = false;
      for (long unsigned int jj = 0; jj < NewList.size(); jj++) {
        if (OldList[ii] == NewList[jj]) {
          RepeatEntry = true;
        }
      }
      if (!RepeatEntry) {
        NewList.push_back(OldList[ii]);
      }
    }
  }

  template <typename T>
  bool checkTrackSelection(T track, std::vector<int64_t> DauParticleList, long unsigned int& skippingPosition, int& rejectionTag)
  {
    if (track.tpcNClsCrossedRows() < 70) {
      rejectionTag = 1;
      return false;
    }
    if (fabs(track.dcaXY()) > 0.2) {
      rejectionTag = 2;
      return false;
    }
    if (!track.isGlobalTrack()) {
      rejectionTag = 3;
      return false;
    }

    bool FlagDaughterTrack = false;
    if (track.globalIndex() == DauParticleList[skippingPosition]) {
      FlagDaughterTrack = true;
      skippingPosition++;
      // LOG(info) <<"*******************iEvent = "<<iEvent<<" :: ID = "<<collIdx<<" :: index skipping = "<<track.globalIndex();
    }
    if (FlagDaughterTrack) {
      rejectionTag = 4;
      return false;
    }
    // Event and Pt filter will filter out some tracks and collisions used for v0 reconstruction.
    if (track.globalIndex() > DauParticleList[skippingPosition] && skippingPosition < DauParticleList.size()) {
      skippingPosition++;
    }
    return true;
  }

  template <typename T>
  void FillFullTrackQA(T track)
  {
    // FullTrack Information
    recoTracks.fill(HIST("htracks_00_1_0_FullTrack_P"), track.p());
    recoTracks.fill(HIST("htracks_00_1_0_FullTrack_tpcInnerParam"), track.tpcInnerParam());
    recoTracks.fill(HIST("htracks_00_1_0_FullTrack_tofExpMom"), track.tofExpMom());

    recoTracks.fill(HIST("hTracks_00_1_1_FullTrack_Pt"), track.pt());
    recoTracks.fill(HIST("hTracks_00_1_2_FullTrack_Eta"), track.eta());
    recoTracks.fill(HIST("hTracks_00_1_3_FullTrack_Phi"), track.phi());
    recoTracks.fill(HIST("hTracks_00_1_4_FullTrack_DcaXY"), track.dcaXY());
    recoTracks.fill(HIST("hTracks_00_1_5_FullTrack_DcaZ"), track.dcaZ());
    recoTracks.fill(HIST("hTracks_00_1_6_FullTrack_Sign"), track.sign());

    // DcaXY
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_00_DcaXY"), track.p(), track.dcaXY());
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_01_DcaXY"), track.pt(), track.dcaXY());
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_02_DcaXY"), track.tpcInnerParam(), track.dcaXY());
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_03_DcaXY"), track.tofExpMom(), track.dcaXY());

    // DcaZ
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_00_DcaZ"), track.p(), track.dcaZ());
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_01_DcaZ"), track.pt(), track.dcaZ());
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_02_DcaZ"), track.tpcInnerParam(), track.dcaZ());
    recoTracks.fill(HIST("htracks_00_1_7_FullTrack_03_DcaZ"), track.tofExpMom(), track.dcaZ());

    // momemtum
    recoTracks.fill(HIST("hTracks_00_2_1_FullTrack_Axis_00_01"), track.p(), track.pt());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_2_2_FullTrack_Axis_00_02"), track.p(), track.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_2_3_FullTrack_Axis_00_03"), track.p(), track.tofExpMom());     // Axis_tofExpMom     ;

    // tpcSignal
    recoTracks.fill(HIST("hTracks_00_3_1_FullTrack_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_3_2_FullTrack_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_3_3_FullTrack_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;

    // tofBeta
    recoTracks.fill(HIST("hTracks_00_4_1_FullTrack_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_4_2_FullTrack_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_4_3_FullTrack_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

    // Look at Pion
    recoTracks.fill(HIST("hTracks_00_5_1_FullTrack_Axis_00_20"), track.p(), track.tpcNSigmaPi());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_5_2_FullTrack_Axis_01_20"), track.pt(), track.tpcNSigmaPi());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_5_3_FullTrack_Axis_02_20"), track.tpcInnerParam(), track.tpcNSigmaPi()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_5_4_FullTrack_Axis_03_20"), track.tofExpMom(), track.tpcNSigmaPi());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_5_5_FullTrack_Axis_00_21"), track.p(), track.tofNSigmaPi());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_5_6_FullTrack_Axis_01_21"), track.pt(), track.tofNSigmaPi());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_5_7_FullTrack_Axis_02_21"), track.tpcInnerParam(), track.tofNSigmaPi()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_5_8_FullTrack_Axis_03_21"), track.tofExpMom(), track.tofNSigmaPi());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_5_9_FullTrack_Axis_20_21"), track.tpcNSigmaPi(), track.tofNSigmaPi());   // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    recoTracks.fill(HIST("hTracks_00_6_1_FullTrack_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_6_2_FullTrack_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_6_3_FullTrack_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_6_4_FullTrack_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_6_5_FullTrack_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_6_6_FullTrack_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_6_7_FullTrack_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_6_8_FullTrack_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_6_9_FullTrack_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    recoTracks.fill(HIST("hTracks_00_7_1_FullTrack_Axis_00_24"), track.p(), track.tpcNSigmaPr());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_7_2_FullTrack_Axis_01_24"), track.pt(), track.tpcNSigmaPr());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_7_3_FullTrack_Axis_02_24"), track.tpcInnerParam(), track.tpcNSigmaPr()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_7_4_FullTrack_Axis_03_24"), track.tofExpMom(), track.tpcNSigmaPr());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_7_5_FullTrack_Axis_00_25"), track.p(), track.tofNSigmaPr());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_7_6_FullTrack_Axis_01_25"), track.pt(), track.tofNSigmaPr());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_7_7_FullTrack_Axis_02_25"), track.tpcInnerParam(), track.tofNSigmaPr()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_7_8_FullTrack_Axis_03_25"), track.tofExpMom(), track.tofNSigmaPr());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_7_9_FullTrack_Axis_24_25"), track.tpcNSigmaPr(), track.tofNSigmaPr());   // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    recoTracks.fill(HIST("hTracks_00_8_1_FullTrack_Axis_00_26"), track.p(), track.tpcNSigmaEl());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_8_2_FullTrack_Axis_01_26"), track.pt(), track.tpcNSigmaEl());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_8_3_FullTrack_Axis_02_26"), track.tpcInnerParam(), track.tpcNSigmaEl()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_8_4_FullTrack_Axis_03_26"), track.tofExpMom(), track.tpcNSigmaEl());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_8_5_FullTrack_Axis_00_27"), track.p(), track.tofNSigmaEl());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_8_6_FullTrack_Axis_01_27"), track.pt(), track.tofNSigmaEl());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_8_7_FullTrack_Axis_02_27"), track.tpcInnerParam(), track.tofNSigmaEl()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_8_8_FullTrack_Axis_03_27"), track.tofExpMom(), track.tofNSigmaEl());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_8_9_FullTrack_Axis_26_27"), track.tpcNSigmaEl(), track.tofNSigmaEl());   // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    recoTracks.fill(HIST("hTracks_00_9_1_FullTrack_Axis_00_28"), track.p(), track.tpcNSigmaDe());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_9_2_FullTrack_Axis_01_28"), track.pt(), track.tpcNSigmaDe());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_9_3_FullTrack_Axis_02_28"), track.tpcInnerParam(), track.tpcNSigmaDe()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_9_4_FullTrack_Axis_03_28"), track.tofExpMom(), track.tpcNSigmaDe());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_9_5_FullTrack_Axis_00_29"), track.p(), track.tofNSigmaDe());             // Axis_p             ;
    recoTracks.fill(HIST("hTracks_00_9_6_FullTrack_Axis_01_29"), track.pt(), track.tofNSigmaDe());            // Axis_pt            ;
    recoTracks.fill(HIST("hTracks_00_9_7_FullTrack_Axis_02_29"), track.tpcInnerParam(), track.tofNSigmaDe()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("hTracks_00_9_8_FullTrack_Axis_03_29"), track.tofExpMom(), track.tofNSigmaDe());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("hTracks_00_9_9_FullTrack_Axis_28_29"), track.tpcNSigmaDe(), track.tofNSigmaDe());   // Axis_tpcInnerParam ;
    // Deuteron
    //
  }

  template <typename T>
  void FillSelectedTrackQA(T track)
  {
    // Full Track Information
    recoTracks.fill(HIST("htracks_11_1_0_SelectedTrack_P"), track.p());
    recoTracks.fill(HIST("htracks_11_1_0_SelectedTrack_tpcInnerParam"), track.tpcInnerParam());
    recoTracks.fill(HIST("htracks_11_1_0_SelectedTrack_tofExpMom"), track.tofExpMom());

    recoTracks.fill(HIST("htracks_11_1_1_SelectedTrack_Pt"), track.pt());
    recoTracks.fill(HIST("htracks_11_1_2_SelectedTrack_Eta"), track.eta());
    recoTracks.fill(HIST("htracks_11_1_3_SelectedTrack_Phi"), track.phi());
    recoTracks.fill(HIST("htracks_11_1_4_SelectedTrack_DcaXY"), track.dcaXY());
    recoTracks.fill(HIST("htracks_11_1_5_SelectedTrack_DcaZ"), track.dcaZ());
    recoTracks.fill(HIST("htracks_11_1_6_SelectedTrack_Sign"), track.sign());

    // DcaXY
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_00_DcaXY"), track.p(), track.dcaXY());
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_01_DcaXY"), track.pt(), track.dcaXY());
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_02_DcaXY"), track.tpcInnerParam(), track.dcaXY());
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_03_DcaXY"), track.tofExpMom(), track.dcaXY());

    // DcaZ
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_00_DcaZ"), track.p(), track.dcaZ());
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_01_DcaZ"), track.pt(), track.dcaZ());
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_02_DcaZ"), track.tpcInnerParam(), track.dcaZ());
    recoTracks.fill(HIST("htracks_11_1_7_SelectedTrack_03_DcaZ"), track.tofExpMom(), track.dcaZ());

    // momemtum
    recoTracks.fill(HIST("htracks_11_2_1_SelectedTrack_Axis_00_01"), track.p(), track.pt());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_2_2_SelectedTrack_Axis_00_02"), track.p(), track.tpcInnerParam()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_2_3_SelectedTrack_Axis_00_03"), track.p(), track.tofExpMom());     // Axis_tofExpMom     ;

    // tpcSignal
    recoTracks.fill(HIST("htracks_11_3_1_SelectedTrack_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_3_2_SelectedTrack_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_3_3_SelectedTrack_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;

    // tofBeta
    recoTracks.fill(HIST("htracks_11_4_1_SelectedTrack_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_4_2_SelectedTrack_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_4_3_SelectedTrack_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

    // Look at Pion
    recoTracks.fill(HIST("htracks_11_5_1_SelectedTrack_Axis_00_20"), track.p(), track.tpcNSigmaPi());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_5_2_SelectedTrack_Axis_01_20"), track.pt(), track.tpcNSigmaPi());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_5_3_SelectedTrack_Axis_02_20"), track.tpcInnerParam(), track.tpcNSigmaPi()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_5_4_SelectedTrack_Axis_03_20"), track.tofExpMom(), track.tpcNSigmaPi());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_5_5_SelectedTrack_Axis_00_21"), track.p(), track.tofNSigmaPi());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_5_6_SelectedTrack_Axis_01_21"), track.pt(), track.tofNSigmaPi());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_5_7_SelectedTrack_Axis_02_21"), track.tpcInnerParam(), track.tofNSigmaPi()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_5_8_SelectedTrack_Axis_03_21"), track.tofExpMom(), track.tofNSigmaPi());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_5_9_SelectedTrack_Axis_20_21"), track.tpcNSigmaPi(), track.tofNSigmaPi());   // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    recoTracks.fill(HIST("htracks_11_6_1_SelectedTrack_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_6_2_SelectedTrack_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_6_3_SelectedTrack_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_6_4_SelectedTrack_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_6_5_SelectedTrack_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_6_6_SelectedTrack_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_6_7_SelectedTrack_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_6_8_SelectedTrack_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_6_9_SelectedTrack_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    recoTracks.fill(HIST("htracks_11_7_1_SelectedTrack_Axis_00_24"), track.p(), track.tpcNSigmaPr());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_7_2_SelectedTrack_Axis_01_24"), track.pt(), track.tpcNSigmaPr());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_7_3_SelectedTrack_Axis_02_24"), track.tpcInnerParam(), track.tpcNSigmaPr()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_7_4_SelectedTrack_Axis_03_24"), track.tofExpMom(), track.tpcNSigmaPr());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_7_5_SelectedTrack_Axis_00_25"), track.p(), track.tofNSigmaPr());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_7_6_SelectedTrack_Axis_01_25"), track.pt(), track.tofNSigmaPr());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_7_7_SelectedTrack_Axis_02_25"), track.tpcInnerParam(), track.tofNSigmaPr()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_7_8_SelectedTrack_Axis_03_25"), track.tofExpMom(), track.tofNSigmaPr());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_7_9_SelectedTrack_Axis_24_25"), track.tpcNSigmaPr(), track.tofNSigmaPr());   // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    recoTracks.fill(HIST("htracks_11_8_1_SelectedTrack_Axis_00_26"), track.p(), track.tpcNSigmaEl());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_8_2_SelectedTrack_Axis_01_26"), track.pt(), track.tpcNSigmaEl());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_8_3_SelectedTrack_Axis_02_26"), track.tpcInnerParam(), track.tpcNSigmaEl()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_8_4_SelectedTrack_Axis_03_26"), track.tofExpMom(), track.tpcNSigmaEl());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_8_5_SelectedTrack_Axis_00_27"), track.p(), track.tofNSigmaEl());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_8_6_SelectedTrack_Axis_01_27"), track.pt(), track.tofNSigmaEl());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_8_7_SelectedTrack_Axis_02_27"), track.tpcInnerParam(), track.tofNSigmaEl()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_8_8_SelectedTrack_Axis_03_27"), track.tofExpMom(), track.tofNSigmaEl());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_8_9_SelectedTrack_Axis_26_27"), track.tpcNSigmaEl(), track.tofNSigmaEl());   // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    recoTracks.fill(HIST("htracks_11_9_1_SelectedTrack_Axis_00_28"), track.p(), track.tpcNSigmaDe());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_9_2_SelectedTrack_Axis_01_28"), track.pt(), track.tpcNSigmaDe());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_9_3_SelectedTrack_Axis_02_28"), track.tpcInnerParam(), track.tpcNSigmaDe()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_9_4_SelectedTrack_Axis_03_28"), track.tofExpMom(), track.tpcNSigmaDe());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_9_5_SelectedTrack_Axis_00_29"), track.p(), track.tofNSigmaDe());             // Axis_p             ;
    recoTracks.fill(HIST("htracks_11_9_6_SelectedTrack_Axis_01_29"), track.pt(), track.tofNSigmaDe());            // Axis_pt            ;
    recoTracks.fill(HIST("htracks_11_9_7_SelectedTrack_Axis_02_29"), track.tpcInnerParam(), track.tofNSigmaDe()); // Axis_tpcInnerParam ;
    recoTracks.fill(HIST("htracks_11_9_8_SelectedTrack_Axis_03_29"), track.tofExpMom(), track.tofNSigmaDe());     // Axis_tofExpMom     ;
    recoTracks.fill(HIST("htracks_11_9_9_SelectedTrack_Axis_28_29"), track.tpcNSigmaDe(), track.tofNSigmaDe());   // Axis_tpcInnerParam ;
    // Deuteron
    //
  }

  template <typename T>
  void fillPionQA(T track, int IdMethod)
  {
    if (IdMethod == 0) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_1_Axis_00_20"), track.p(), track.tpcNSigmaPi());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_2_Axis_01_20"), track.pt(), track.tpcNSigmaPi());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_3_Axis_02_20"), track.tpcInnerParam(), track.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_4_Axis_03_20"), track.tofExpMom(), track.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_5_Axis_00_21"), track.p(), track.tofNSigmaPi());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_6_Axis_01_21"), track.pt(), track.tofNSigmaPi());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_7_Axis_02_21"), track.tpcInnerParam(), track.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_8_Axis_03_21"), track.tofExpMom(), track.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id0_3_9_Axis_20_21"), track.tpcNSigmaPi(), track.tofNSigmaPi());   // Axis_tpcInnerParam ;
    } else if (IdMethod == 1) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_1_Axis_00_20"), track.p(), track.tpcNSigmaPi());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_2_Axis_01_20"), track.pt(), track.tpcNSigmaPi());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_3_Axis_02_20"), track.tpcInnerParam(), track.tpcNSigmaPi()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_4_Axis_03_20"), track.tofExpMom(), track.tpcNSigmaPi());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_5_Axis_00_21"), track.p(), track.tofNSigmaPi());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_6_Axis_01_21"), track.pt(), track.tofNSigmaPi());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_7_Axis_02_21"), track.tpcInnerParam(), track.tofNSigmaPi()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_8_Axis_03_21"), track.tofExpMom(), track.tofNSigmaPi());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_11_Pi_Id1_3_9_Axis_20_21"), track.tpcNSigmaPi(), track.tofNSigmaPi());   // Axis_tpcInnerParam ;
    }
  }

  template <typename T>
  void fillKaonQA(T track, int IdMethod)
  {
    if (IdMethod == 0) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_1_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_2_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_3_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_4_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_5_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_6_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_7_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_8_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id0_3_9_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    } else if (IdMethod == 1) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_1_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_2_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_3_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_4_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_5_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_6_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_7_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_8_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_12_Ka_Id1_3_9_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    }
  }

  template <typename T>
  void fillProtonQA(T track, int IdMethod)
  {
    if (IdMethod == 0) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_1_Axis_00_24"), track.p(), track.tpcNSigmaPr());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_2_Axis_01_24"), track.pt(), track.tpcNSigmaPr());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_3_Axis_02_24"), track.tpcInnerParam(), track.tpcNSigmaPr()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_4_Axis_03_24"), track.tofExpMom(), track.tpcNSigmaPr());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_5_Axis_00_25"), track.p(), track.tofNSigmaPr());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_6_Axis_01_25"), track.pt(), track.tofNSigmaPr());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_7_Axis_02_25"), track.tpcInnerParam(), track.tofNSigmaPr()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_8_Axis_03_25"), track.tofExpMom(), track.tofNSigmaPr());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id0_3_9_Axis_24_25"), track.tpcNSigmaPr(), track.tofNSigmaPr());   // Axis_tpcInnerParam ;
    } else if (IdMethod == 1) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_1_Axis_00_24"), track.p(), track.tpcNSigmaPr());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_2_Axis_01_24"), track.pt(), track.tpcNSigmaPr());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_3_Axis_02_24"), track.tpcInnerParam(), track.tpcNSigmaPr()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_4_Axis_03_24"), track.tofExpMom(), track.tpcNSigmaPr());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_5_Axis_00_25"), track.p(), track.tofNSigmaPr());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_6_Axis_01_25"), track.pt(), track.tofNSigmaPr());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_7_Axis_02_25"), track.tpcInnerParam(), track.tofNSigmaPr()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_8_Axis_03_25"), track.tofExpMom(), track.tofNSigmaPr());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_13_Pr_Id1_3_9_Axis_24_25"), track.tpcNSigmaPr(), track.tofNSigmaPr());   // Axis_tpcInnerParam ;
    }
  }

  template <typename T>
  void fillElectronQA(T track, int IdMethod)
  {
    if (IdMethod == 0) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_1_Axis_00_26"), track.p(), track.tpcNSigmaEl());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_2_Axis_01_26"), track.pt(), track.tpcNSigmaEl());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_3_Axis_02_26"), track.tpcInnerParam(), track.tpcNSigmaEl()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_4_Axis_03_26"), track.tofExpMom(), track.tpcNSigmaEl());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_5_Axis_00_27"), track.p(), track.tofNSigmaEl());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_6_Axis_01_27"), track.pt(), track.tofNSigmaEl());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_7_Axis_02_27"), track.tpcInnerParam(), track.tofNSigmaEl()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_8_Axis_03_27"), track.tofExpMom(), track.tofNSigmaEl());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id0_3_9_Axis_26_27"), track.tpcNSigmaEl(), track.tofNSigmaEl());   // Axis_tpcInnerParam ;
    } else if (IdMethod == 1) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_1_Axis_00_26"), track.p(), track.tpcNSigmaEl());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_2_Axis_01_26"), track.pt(), track.tpcNSigmaEl());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_3_Axis_02_26"), track.tpcInnerParam(), track.tpcNSigmaEl()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_4_Axis_03_26"), track.tofExpMom(), track.tpcNSigmaEl());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_5_Axis_00_27"), track.p(), track.tofNSigmaEl());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_6_Axis_01_27"), track.pt(), track.tofNSigmaEl());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_7_Axis_02_27"), track.tpcInnerParam(), track.tofNSigmaEl()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_8_Axis_03_27"), track.tofExpMom(), track.tofNSigmaEl());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_14_El_Id1_3_9_Axis_26_27"), track.tpcNSigmaEl(), track.tofNSigmaEl());   // Axis_tpcInnerParam ;
    }
  }

  template <typename T>
  void fillDeuteronQA(T track, int IdMethod)
  {
    if (IdMethod == 0) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_1_Axis_00_28"), track.p(), track.tpcNSigmaDe());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_2_Axis_01_28"), track.pt(), track.tpcNSigmaDe());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_3_Axis_02_28"), track.tpcInnerParam(), track.tpcNSigmaDe()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_4_Axis_03_28"), track.tofExpMom(), track.tpcNSigmaDe());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_5_Axis_00_29"), track.p(), track.tofNSigmaDe());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_6_Axis_01_29"), track.pt(), track.tofNSigmaDe());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_7_Axis_02_29"), track.tpcInnerParam(), track.tofNSigmaDe()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_8_Axis_03_29"), track.tofExpMom(), track.tofNSigmaDe());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id0_3_9_Axis_28_29"), track.tpcNSigmaDe(), track.tofNSigmaDe());   // Axis_tpcInnerParam ;
    } else if (IdMethod == 1) {
      // tpcSignal
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_1_Axis_00_28"), track.p(), track.tpcNSigmaDe());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_2_Axis_01_28"), track.pt(), track.tpcNSigmaDe());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_3_Axis_02_28"), track.tpcInnerParam(), track.tpcNSigmaDe()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_4_Axis_03_28"), track.tofExpMom(), track.tpcNSigmaDe());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_5_Axis_00_29"), track.p(), track.tofNSigmaDe());             // Axis_p             ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_6_Axis_01_29"), track.pt(), track.tofNSigmaDe());            // Axis_pt            ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_7_Axis_02_29"), track.tpcInnerParam(), track.tofNSigmaDe()); // Axis_tpcInnerParam ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_8_Axis_03_29"), track.tofExpMom(), track.tofNSigmaDe());     // Axis_tofExpMom     ;
      recoAnalysis.fill(HIST("hAnalysis_15_De_Id1_3_9_Axis_28_29"), track.tpcNSigmaDe(), track.tofNSigmaDe());   // Axis_tpcInnerParam ;
    }
  }

  // Event Filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZvertex);

  // Track Filter
  Filter PtFilter = (o2::aod::track::pt) > 0.15f && (o2::aod::track::pt) < 2.0f;

  // Filters on V0s
  Filter preFilterv0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                               aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>; // aod::CentFV0As,

  using myTracks = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>>;

  // For manual sliceBy
  Preslice<myTracks> TracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<aod::V0Datas> V0sPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<aod::CascDatas> CascadesPerCollisionPreslice = o2::aod::track::collisionId;

  int CountDebug = 0;
  int checker = 0;
  int iEvent = -1;
  double lowest = 0;
  double maximum = 0;
  int dfcount = 0;
  bool UnusedError = false;
  int v0count = 0;
  void process(myCollisions const& collisions,
               soa::Filtered<aod::V0Datas> const& V0s,
               myTracks const& tracks) // aod::TracksIU const& FullTrack)
  {

    std::vector<int64_t> K0sPosList; // long int
    std::vector<int64_t> K0sNegList; // long int

    std::vector<std::vector<int64_t>> K0sPosListWithCollision; // long int
    std::vector<std::vector<int64_t>> K0sNegListWithCollision; // long int

    // iterate over all the v0s of current dataframe
    for (const auto& v0 : V0s) {
      recoV0s.fill(HIST("hV0s_0_0_01_K0s_Mass_Full"), v0.mK0Short());
      recoV0s.fill(HIST("hV0s_0_0_02_Lambda_Mass_Full"), v0.mLambda());
      recoV0s.fill(HIST("hV0s_0_0_03_AntiLambda_Mass_Full"), v0.mAntiLambda());

      // cut on dynamic columns for v0 particles
      if (v0.v0cosPA() < v0setting_cospa)
        continue; // collision.posX(), collision.posY(), collision.posZ()
      if (v0.v0radius() < v0setting_radius)
        continue;

      const auto& posDaughterTrack = v0.posTrack_as<myTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<myTracks>();

      int posIdMethod = -1;
      int negIdMethod = -1;

      int posPiIdMethod = -1;
      int posPrIdMethod = -1;
      int negPiIdMethod = -1;
      int negPrIdMethod = -1;
      int V0Tag = FindV0Tag(posDaughterTrack, negDaughterTrack, posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod);
      // v0count++;
      // LOG(info)<<"v0count = "<<v0count<<" ::OUTER:: posIdMethod = "<<posIdMethod<<" :: negIdMethod = "<<negIdMethod;
      bool isK0s = false;
      bool isLambda = false;
      bool isAntiLambda = false;

      if ((V0Tag & 1) == 1)
        isK0s = true;
      if ((V0Tag & 2) == 2)
        isLambda = true;
      if ((V0Tag & 4) == 4)
        isAntiLambda = true;

      if (isLambda == isAntiLambda && isLambda != isAntiLambda) {
        LOG(info) << " Do Nothing : Just Using the Variables to bypass the errors";
      }

      int TrueV0TagValue = 0;
      // K0s Analysis
      if (isK0s) {
        posIdMethod = posPiIdMethod;
        negIdMethod = negPiIdMethod;
        // topological information.
        recoV0s.fill(HIST("hV0s_0_1_01_K0s_dcapostopv"), v0.dcapostopv());
        recoV0s.fill(HIST("hV0s_0_1_02_K0s_dcanegtopv"), v0.dcanegtopv());
        recoV0s.fill(HIST("hV0s_0_1_03_K0s_dcaV0daughters"), v0.dcaV0daughters());
        recoV0s.fill(HIST("hV0s_0_1_04_K0s_v0cosPA"), v0.v0cosPA());
        recoV0s.fill(HIST("hV0s_0_1_05_K0s_v0radius"), v0.v0radius());

        recoV0s.fill(HIST("hV0s_0_2_01_K0s_Mass_DynamicCuts"), v0.mK0Short());

        // K0s only using mass cuts, no pt and eta cuts
        if (mLowK0s < v0.mK0Short() && v0.mK0Short() < mHighK0s) {
          fillV0_K0sFullInformation(v0, posDaughterTrack, negDaughterTrack, posIdMethod, negIdMethod);
        }
        // Select Your Kaon Here.
        if (selK0s(v0)) {
          fillV0_K0sSelectedInformation(v0, posDaughterTrack, negDaughterTrack, posIdMethod, negIdMethod);

          TrueV0TagValue += 1;                                  // nK0s++;
          K0sPosList.push_back(posDaughterTrack.globalIndex()); // posDaughterTrack.index());
          K0sNegList.push_back(negDaughterTrack.globalIndex()); // negDaughterTrack.index());
          K0sPosListWithCollision.push_back({posDaughterTrack.globalIndex(), posDaughterTrack.collisionId()});
          K0sNegListWithCollision.push_back({negDaughterTrack.globalIndex(), negDaughterTrack.collisionId()});
        }
      } // End of K0s loop

      recoV0s.fill(HIST("hV0s_0_0_00_hV0TagCount"), V0Tag);
      recoV0s.fill(HIST("hV0s_0_0_00_hTrueV0TagCount"), TrueV0TagValue);
    } // End of V0s Loop

    // LOG(info)<<"v0 size = "<<V0s.size();
    FindRepeatEntries(K0sPosList, recoV0s.get<TH1>(HIST("hV0s_0_4_02_K0s_nCommonPionOfDifferentK0s")));
    FindRepeatEntries(K0sNegList, recoV0s.get<TH1>(HIST("hV0s_0_4_02_K0s_nCommonPionOfDifferentK0s")));

    std::vector<int64_t> V0PosList;
    std::vector<int64_t> V0NegList;

    FillNewListFromOldList(V0PosList, K0sPosList);
    FillNewListFromOldList(V0NegList, K0sNegList);

    // //Final Track Vectors
    std::vector<int64_t> PosDauList;
    std::vector<int64_t> NegDauList;

    FillNewListFromOldList(PosDauList, V0PosList);
    FillNewListFromOldList(NegDauList, V0NegList);

    // FullDauList
    std::vector<int64_t> FullDauList;
    FillNewListFromOldList(FullDauList, PosDauList);
    FillNewListFromOldList(FullDauList, NegDauList);

    InsertionSortVector(PosDauList);
    InsertionSortVector(FullDauList);

    for (long unsigned int ii = 0; ii < (PosDauList.size() - 1); ii++) {
      if (PosDauList[ii] > PosDauList[ii + 1]) {
        LOG(info) << "ERROR:: PosDauList Not in Sequence = " << ii;
      }
    }

    for (long unsigned int ii = 0; ii < (FullDauList.size() - 1); ii++) {
      if (FullDauList[ii] > FullDauList[ii + 1]) {
        LOG(info) << "ERROR:: FullDauList Not in Sequence = " << ii;
      }
    }

    long unsigned int skippingPosition = 0;
    for (const auto& collision : collisions) {
      iEvent++;
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
      int nPion = 0;
      int nKaon = 0;

      double CentFT0C = collision.centFT0C();
      recoEvent.fill(HIST("hCentrality"), CentFT0C);
      recoEvent.fill(HIST("hCollisionCount"), 0.5);
      recoEvent.fill(HIST("hVertexXRec"), collision.posX());
      recoEvent.fill(HIST("hVertexYRec"), collision.posY());
      recoEvent.fill(HIST("hVertexZRec"), collision.posZ());

      // group tracks, v0s, Cascades manually
      const uint64_t collIdx = collision.globalIndex();
      const auto TracksTable_perColl = tracks.sliceBy(TracksPerCollisionPreslice, collIdx);
      const auto V0sTable_perColl = V0s.sliceBy(V0sPerCollisionPreslice, collIdx);

      recoEvent.fill(HIST("hV0Size"), V0sTable_perColl.size());

      nK0s = 0;
      // centrality dependent mass
      for (const auto& v0 : V0sTable_perColl) {
        recoK0s.fill(HIST("hK0s_2_00_Mass_Full"), v0.mK0Short());
        recoLambda.fill(HIST("hLambda_Mass_Full"), v0.mLambda());
        recoLambda.fill(HIST("hAntiLambda_Mass_Full"), v0.mAntiLambda());

        // cut on dynamic columns for v0 particles
        if (v0.v0cosPA() < v0setting_cospa)
          continue; // collision.posX(), collision.posY(), collision.posZ()
        if (v0.v0radius() < v0setting_radius)
          continue;

        const auto& posDaughterTrack = v0.posTrack_as<myTracks>();
        const auto& negDaughterTrack = v0.negTrack_as<myTracks>();

        int posIdMethod = -1;
        int negIdMethod = -1;

        bool isK0s = false;
        bool isLambda = false;
        bool isAntiLambda = false;

        int posPiIdMethod = -1;
        int posPrIdMethod = -1;
        int negPiIdMethod = -1;
        int negPrIdMethod = -1;
        int V0Tag = FindV0Tag(posDaughterTrack, negDaughterTrack, posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod);

        if ((V0Tag & 1) == 1)
          isK0s = true;
        if ((V0Tag & 2) == 2)
          isLambda = true;
        if ((V0Tag & 4) == 4)
          isAntiLambda = true;

        if (isLambda == isAntiLambda && isLambda != isAntiLambda) {
          LOG(info) << " Do Nothing : Just Using the Variables to bypass the errors";
        }

        // K0s Analysis
        if (isK0s) {
          posIdMethod = posPiIdMethod;
          negIdMethod = negPiIdMethod;
          // topological information.
          recoK0s.fill(HIST("hK0s_1_01_dcapostopv"), v0.dcapostopv());
          recoK0s.fill(HIST("hK0s_1_02_dcanegtopv"), v0.dcanegtopv());
          recoK0s.fill(HIST("hK0s_1_03_dcaV0daughters"), v0.dcaV0daughters());
          recoK0s.fill(HIST("hK0s_1_04_v0cosPA"), v0.v0cosPA());
          recoK0s.fill(HIST("hK0s_1_05_v0radius"), v0.v0radius());

          recoK0s.fill(HIST("hK0s_2_01_Mass_DynamicCuts"), v0.mK0Short());

          if (selK0s(v0)) {
            fillK0sSelectedInformation(v0, posDaughterTrack, negDaughterTrack, posIdMethod, negIdMethod);
            recoK0s.fill(HIST("hK0s_5_01_centFTOC_mK0s"), CentFT0C, v0.mK0Short());
            nK0s++;
          }
        } // End of K0s loop
      } // End of V0s Loop

      nTrack = 0;
      int nRejectedPiPlus = 0;
      int nRejectedPiMinus = 0;
      for (const auto& track : TracksTable_perColl) {

        FillFullTrackQA(track);
        int rejectionTag = 0;
        if (!checkTrackSelection(track, FullDauList, skippingPosition, rejectionTag)) {
          if (rejectionTag == 4) {
            if (track.sign() > 0) {
              nRejectedPiPlus++;
            }
            if (track.sign() < 0) {
              nRejectedPiMinus++;
            }
          }
          recoAnalysis.fill(HIST("hAnalysis_4_RejectedTrack_RejectionTag"), rejectionTag);
          continue;
        }

        FillSelectedTrackQA(track);

        nTrack++;
        // Do Proper Track Identification
        bool TrackIsPion = false;
        bool TrackIsKaon = false;
        bool TrackIsProton = false;
        bool TrackIsElectron = false;
        bool TrackIsDeuteron = false;

        int TrackIdTag = 0;
        int PiIdMethod = -1;
        int KaIdMethod = -1;
        int PrIdMethod = -1;
        int ElIdMethod = -1;
        int DeIdMethod = -1;

        if (selPion(track, PiIdMethod)) {
          TrackIsPion = true;
          TrackIdTag += 1;
        }
        if (selKaon(track, KaIdMethod)) {
          TrackIsKaon = true;
          TrackIdTag += 2;
        }
        if (selProton(track, PrIdMethod)) {
          TrackIsProton = true;
          TrackIdTag += 4;
        }
        if (selElectron(track, ElIdMethod)) {
          TrackIsElectron = true;
          TrackIdTag += 8;
        }
        if (selDeuteron(track, DeIdMethod)) {
          TrackIsDeuteron = true;
          TrackIdTag += 16;
        }

        if (TrackIsPion) {
          fillPionQA(track, PiIdMethod);
          if (track.sign() > 0) {
            nPiPlus++;
          }
          if (track.sign() < 0) {
            nPiMinus++;
          }
        }
        if (TrackIsKaon) {
          fillKaonQA(track, KaIdMethod);
          if (track.sign() > 0) {
            nKaPlus++;
          }
          if (track.sign() < 0) {
            nKaMinus++;
          }
        }
        if (TrackIsProton) {
          fillProtonQA(track, PrIdMethod);
          if (track.sign() > 0) {
            nProton++;
          }
          if (track.sign() < 0) {
            nPBar++;
          }
        }
        if (TrackIsElectron) {
          fillElectronQA(track, ElIdMethod);
          if (track.sign() > 0) {
            nElPlus++;
          }
          if (track.sign() < 0) {
            nElMinus++;
          }
        }
        if (TrackIsDeuteron) {
          fillDeuteronQA(track, DeIdMethod);
          if (track.sign() > 0) {
            nDePlus++;
          }
          if (track.sign() < 0) {
            nDeMinus++;
          }
        }
        // hAnalysis_4_SelectedTrack_IdentificationTag
        recoAnalysis.fill(HIST("hAnalysis_4_SelectedTrack_IdentificationTag"), TrackIdTag);
      } // track loop ends

      nPion = nPiPlus + nPiMinus;
      nKaon = nKaPlus + nKaMinus;

      recoAnalysis.fill(HIST("hAnalysis_5_Sparse_Full_K0sPiKa"),
                        CentFT0C,
                        nTrack,
                        nK0s,
                        nRejectedPiPlus,
                        nRejectedPiMinus,
                        nPiPlus,
                        nPiMinus,
                        nKaPlus,
                        nKaMinus);
      recoAnalysis.fill(HIST("hAnalysis_5_Sparse_Full_K0sPrDe"),
                        CentFT0C,
                        nTrack,
                        nK0s,
                        nRejectedPiPlus,
                        nRejectedPiMinus,
                        nProton,
                        nPBar,
                        nDePlus,
                        nDeMinus);
      recoAnalysis.fill(HIST("hAnalysis_5_Sparse_Full_K0sKaEl"),
                        CentFT0C,
                        nTrack,
                        nK0s,
                        nRejectedPiPlus,
                        nRejectedPiMinus,
                        nKaPlus,
                        nKaMinus,
                        nElPlus,
                        nElMinus);
      recoAnalysis.fill(HIST("hAnalysis_5_Sparse_Full_PiKaPr"),
                        CentFT0C,
                        nTrack,
                        nRejectedPiPlus,
                        nRejectedPiMinus,
                        nPiPlus,
                        nPiMinus,
                        nKaPlus,
                        nKaMinus,
                        nProton,
                        nPBar);
      recoAnalysis.fill(HIST("hAnalysis_5_Sparse_Full_PiElDe"),
                        CentFT0C,
                        nTrack,
                        nRejectedPiPlus,
                        nRejectedPiMinus,
                        nPiPlus,
                        nPiMinus,
                        nElPlus,
                        nElMinus,
                        nDePlus,
                        nDeMinus);
      recoAnalysis.fill(HIST("hAnalysis_5_Sparse_Full_KaPrDe"),
                        CentFT0C,
                        nTrack,
                        nRejectedPiPlus,
                        nRejectedPiMinus,
                        nKaPlus,
                        nKaMinus,
                        nProton,
                        nPBar,
                        nDePlus,
                        nDeMinus);
      recoAnalysis.fill(HIST("hAnalysis_5_Sparse_Full_PrElDe"),
                        CentFT0C,
                        nTrack,
                        nRejectedPiPlus,
                        nRejectedPiMinus,
                        nProton,
                        nPBar,
                        nElPlus,
                        nElMinus,
                        nDePlus,
                        nDeMinus);

      if (nK0s > 0) {
        recoAnalysis.fill(HIST("hAnalysis_5_Sparse_newDynm_K0s_Ka"),
                          CentFT0C, nTrack, nK0s, nKaon,
                          nK0s * nK0s, nKaon * nKaon, nK0s * nKaon);
        //
        recoAnalysis.fill(HIST("hAnalysis_5_Sparse_newDynm_Kp_Km"),
                          CentFT0C, nTrack, nKaPlus, nKaMinus,
                          nKaPlus * nKaPlus, nKaMinus * nKaMinus, nKaPlus * nKaMinus);
        //
      }
      recoAnalysis.fill(HIST("hTest"), nPion * nPion);
    } // collision loop ends
  } // Process Function Ends
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<isospin_fluctuation>(cfgc)};
}
