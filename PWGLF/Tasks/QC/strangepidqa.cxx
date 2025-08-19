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
// This task is designed to do QA to the TOF PID applied to strangeness
// in the regular framework

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using std::cout;
using std::endl;

struct strangepidqa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};

  // base properties
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisRadius{"axisRadius", {200, 0.0f, 100.0f}, "V0 radius (cm)"};

  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 50.0f, 70.0f, 100.0f}, "FT0C centrality"};

  AxisSpec massAxisXi = {200, 1.222f, 1.422f, "Inv. Mass (GeV/c^{2})"};
  AxisSpec massAxisOmega = {200, 1.572f, 1.772f, "Inv. Mass (GeV/c^{2})"};

  // Invariant Mass
  ConfigurableAxis axisK0ShortMass{"axisK0ShortMass", {200, 0.497f - 0.050f, 0.497f + 0.050f}, "M_{K0s} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.08f, 1.16f}, "M_{#Lambda} (GeV/c^{2})"};

  // time axes
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {200, -1000.0f, +1000.0f}, "#Delta time (ps)"};
  ConfigurableAxis axisTime{"axisTime", {200, 0.0f, +20000.0f}, "T (ps)"};
  ConfigurableAxis axisBeta{"axisBeta", {1200, 0.0f, +1.2f}, "#Beta"};

  // Length axis
  ConfigurableAxis axisLength{"axisLength", {600, 0.0f, +600.0f}, "track Length (cm)"};

  // TOF cut axis
  ConfigurableAxis axisTOFCut{"axisTOFCut", {100, 0.0f, +10000.0f}, "TOF compat. cut (ps)"};

  // TOF selection matters
  Configurable<bool> requireTOFsignalPion{"requireTOFsignalPion", true, "require that pion prongs have TOF"};
  Configurable<bool> requireTOFsignalProton{"requireTOFsignalProton", true, "require that proton prongs have TOF"};
  Configurable<bool> requireTOFEventTimePion{"requireTOFEventTimePion", true, "require that pion prongs have TOF event time"};
  Configurable<bool> requireTOFEventTimeProton{"requireTOFEventTimeProton", true, "require that proton prongs have TOF event time"};

  Configurable<float> maxDeltaTimeProton{"maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
  Configurable<float> maxDeltaTimePion{"maxDeltaTimePion", 1e+9, "check maximum allowed time"};
  Configurable<float> maxDeltaTimeDecay{"maxDeltaTimeDecay", 1e+9, "check maximum allowed delta-time-decay"};

  Configurable<float> minCentrality{"minCentrality", 60, "max value of centrality allowed"};
  Configurable<float> maxCentrality{"maxCentrality", 100, "max value of centrality allowed"};

  Configurable<float> minV0Radius{"minV0Radius", 1.5, "min radius"};
  Configurable<float> minCosPA{"minCosPA", .98, "min cosPA"};

  Configurable<float> minPtV0{"minPtV0", 1.0, "min pT for integrated mass histograms"};
  Configurable<float> maxPtV0{"maxPtV0", 3.0, "max pT for integrated mass histograms"};

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Selection criteria for cascade analysis
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.95, "v0setting_cospa"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0, "v0setting_dcav0dau"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "v0setting_dcapostopv"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "v0setting_dcanegtopv"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "v0setting_radius"};
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Track configurables
  Configurable<int> tpcCrossedRows{"tpcCrossedRows", 70, "minimum number of TPC rows requirement"};
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 4, "TPC NSigma bachelor (>10 is no cut)"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 4, "TPC NSigma proton <- lambda (>10 is no cut)"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 4, "TPC NSigma pion <- lambda (>10 is no cut)"};

  Configurable<float> tofNsigmaXiLaPr{"tpcNsigmaXiLaPr", 1e+5, "TOF NSigma proton <- lambda <- Xi (>10 is no cut)"};
  Configurable<float> tofNsigmaXiLaPi{"tpcNsigmaXiLaPi", 1e+5, "TOF NSigma pion <- lambda <- Xi (>10 is no cut)"};
  Configurable<float> tofNsigmaXiPi{"tpcNsigmaXiPi", 1e+5, "TOF NSigma pion <- Xi (>10 is no cut)"};
  Configurable<float> tofNsigmaOmLaPr{"tpcNsigmaOmLaPr", 1e+5, "TOF NSigma proton <- lambda <- Omega (>10 is no cut)"};
  Configurable<float> tofNsigmaOmLaPi{"tpcNsigmaOmLaPi", 1e+5, "TOF NSigma pion <- lambda <- Omega (>10 is no cut)"};
  Configurable<float> tofNsigmaOmKa{"tpcNsigmaOmKa", 1e+5, "TOF NSigma Kaon <- Omega (>10 is no cut)"};

  Configurable<float> tofNsigmaCompatibility{"tofNsigmaCompatibility", 4, "compatibility check for V0s"};
  Configurable<float> tofNsigmaCompatibilityCascades{"tofNsigmaCompatibilityCascades", 4, "compatibility check for cascades"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, 100.0f}});

    // Presence of negative and positive signals and event time
    auto hPositiveStatus = histos.add<TH1>("hPositiveStatus", "hPositiveStatus", kTH1F, {{4, -0.5f, 3.5f}});
    auto hNegativeStatus = histos.add<TH1>("hNegativeStatus", "hNegativeStatus", kTH1F, {{4, -0.5f, 3.5f}});
    auto hPositiveStatusReachedTOF = histos.add<TH1>("hPositiveStatusReachedTOF", "hPositiveStatusReachedTOF", kTH1F, {{4, -0.5f, 3.5f}});
    auto hNegativeStatusReachedTOF = histos.add<TH1>("hNegativeStatusReachedTOF", "hNegativeStatusReachedTOF", kTH1F, {{4, -0.5f, 3.5f}});
    hPositiveStatus->GetXaxis()->SetBinLabel(1, "All positive");
    hPositiveStatus->GetXaxis()->SetBinLabel(2, "Has only TOF sig");
    hPositiveStatus->GetXaxis()->SetBinLabel(3, "Has only TOF ev time");
    hPositiveStatus->GetXaxis()->SetBinLabel(4, "Has full time info");
    hNegativeStatus->GetXaxis()->SetBinLabel(1, "All negative");
    hNegativeStatus->GetXaxis()->SetBinLabel(2, "Has only TOF sig");
    hNegativeStatus->GetXaxis()->SetBinLabel(3, "Has only TOF ev time");
    hNegativeStatus->GetXaxis()->SetBinLabel(4, "Has full time info");

    hPositiveStatusReachedTOF->GetXaxis()->SetBinLabel(1, "All positive");
    hPositiveStatusReachedTOF->GetXaxis()->SetBinLabel(2, "Has only TOF sig");
    hPositiveStatusReachedTOF->GetXaxis()->SetBinLabel(3, "Has only TOF ev time");
    hPositiveStatusReachedTOF->GetXaxis()->SetBinLabel(4, "Has full time info");
    hNegativeStatusReachedTOF->GetXaxis()->SetBinLabel(1, "All negative");
    hNegativeStatusReachedTOF->GetXaxis()->SetBinLabel(2, "Has only TOF sig");
    hNegativeStatusReachedTOF->GetXaxis()->SetBinLabel(3, "Has only TOF ev time");
    hNegativeStatusReachedTOF->GetXaxis()->SetBinLabel(4, "Has full time info");

    // V0 Radius
    histos.add("hLambdaMass", "hLambdaMass", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hAssocLambdaMass", "hAssocLambdaMass", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hAssocLambdaMassGoodTime", "hAssocLambdaMassGoodTime", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hAssocLambdaMassBadTime", "hAssocLambdaMassBadTime", {HistType::kTH1F, {axisLambdaMass}});

    // V0 Radius
    histos.add("h2dLambdaRadiusVsPt", "hLambdaRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});

    // Invariant Mass
    histos.add("h2dLambdaMassVsPt", "hLambdaMassVsPt", {HistType::kTH2F, {axisPt, axisLambdaMass}});

    // Invariant Mass
    histos.add("h2dLambdaMassVsTOFCut", "h2dLambdaMassVsTOFCut", {HistType::kTH2F, {axisLambdaMass, axisTOFCut}});
    histos.add("h2dLambdaMassVsTOFCutWithSignals", "h2dLambdaMassVsTOFCutWithSignals", {HistType::kTH2F, {axisLambdaMass, axisTOFCut}});
    histos.add("h2dLambdaMassVsTOFCutMeson", "h2dLambdaMassVsTOFCutMeson", {HistType::kTH2F, {axisLambdaMass, axisTOFCut}});
    histos.add("h2dLambdaMassVsTOFCutMesonWithSignals", "h2dLambdaMassVsTOFCutMesonWithSignals", {HistType::kTH2F, {axisLambdaMass, axisTOFCut}});

    // Invariant Mass with TOF selections
    histos.add("hLambdaMass_ProtonTOF", "hLambdaMass_ProtonTOF", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hLambdaMass_PionTOF", "hLambdaMass_PionTOF", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hLambdaMass_AllTOF", "hLambdaMass_AllTOF", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hLambdaMass_DeltaDecayTime", "hLambdaMass_DeltaDecayTime", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h2dLambdaMassVsPt_ProtonTOF", "hLambdaMassVsPtProtonTOF", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("h2dLambdaMassVsPt_PionTOF", "hLambdaMassVsPtPionTOF", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("h2dLambdaMassVsPt_AllTOF", "hLambdaMassVsPtAllTOF", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("h2dLambdaMassVsPt_DeltaDecayTime", "hLambdaMassVsPtDeltaDecayTime", {HistType::kTH2F, {axisPt, axisLambdaMass}});

    // Invariant Mass with TOF selections
    histos.add("hLambdaMass_InvertProtonTOF", "hLambdaMass_InvertProtonTOF", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hLambdaMass_InvertPionTOF", "hLambdaMass_InvertPionTOF", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("hLambdaMass_InvertAllTOF", "hLambdaMass_InvertAllTOF", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h2dLambdaMassVsPt_InvertProtonTOF", "hLambdaMassVsPtInvertProtonTOF", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("h2dLambdaMassVsPt_InvertPionTOF", "hLambdaMassVsPtInvertPionTOF", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("h2dLambdaMassVsPt_InvertAllTOF", "hLambdaMassVsPtInvertAllTOF", {HistType::kTH2F, {axisPt, axisLambdaMass}});

    // radius vs prong length
    histos.add("h2dProtonLengthVsRadius", "h2dProtonLengthVsRadius", {HistType::kTH2F, {axisRadius, axisLength}});
    histos.add("h2dPionLengthVsRadius", "h2dPionLengthVsRadius", {HistType::kTH2F, {axisRadius, axisLength}});

    // recalculated vs topv lengths
    histos.add("h2dProtonLengthVsLengthToPV", "h2dProtonLengthVsRadiusToPV", {HistType::kTH2F, {axisRadius, axisRadius}});
    histos.add("h2dPionLengthVsLengthToPV", "h2dPionLengthVsLengthToPV", {HistType::kTH2F, {axisRadius, axisRadius}});

    // TOF PID testing for prongs
    histos.add("h2dProtonDeltaTimeVsPt", "h2dProtonDeltaTimeVsPt", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    histos.add("h2dProtonDeltaTimeVsRadius", "h2dProtonDeltaTimeVsRadius", {HistType::kTH2F, {axisRadius, axisDeltaTime}});
    histos.add("h2dPionDeltaTimeVsPt", "h2dPionDeltaTimeVsPt", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    histos.add("h2dPionDeltaTimeVsRadius", "h2dPionDeltaTimeVsRadius", {HistType::kTH2F, {axisRadius, axisDeltaTime}});

    // TOF PID testing for prongs
    histos.add("h2dProtonDeltaTimeVsPt_MassSelected", "h2dProtonDeltaTimeVsPt_MassSelected", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    histos.add("h2dProtonDeltaTimeVsRadius_MassSelected", "h2dProtonDeltaTimeVsRadius_MassSelected", {HistType::kTH2F, {axisRadius, axisDeltaTime}});
    histos.add("h2dPionDeltaTimeVsPt_MassSelected", "h2dPionDeltaTimeVsPt_MassSelected", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    histos.add("h2dPionDeltaTimeVsRadius_MassSelected", "h2dPionDeltaTimeVsRadius_MassSelected", {HistType::kTH2F, {axisRadius, axisDeltaTime}});

    // delta lambda decay time
    histos.add("h2dLambdaDeltaDecayTime", "h2dLambdaDeltaDecayTime", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    histos.add("h2dLambdaDeltaDecayTime_MassSelected", "h2dLambdaDeltaDecayTime_MassSelected", {HistType::kTH2F, {axisPt, axisDeltaTime}});

    // beta plots
    histos.add("h2dLambdaBeta", "h2dLambdaBeta", {HistType::kTH2F, {axisPt, axisBeta}});
    histos.add("h2dLambdaBeta_MassSelected", "h2dLambdaBeta_MassSelected", {HistType::kTH2F, {axisPt, axisBeta}});

    // length vs time for prongs / debug
    histos.add("h2dTimeVsLengthProtonProng", "h2dTimeVsLengthProtonProng", {HistType::kTH2F, {axisLength, axisTime}});
    histos.add("h2dTimeVsLengthPionProng", "h2dTimeVsLengthPionProng", {HistType::kTH2F, {axisLength, axisTime}});

    histos.add("h1dMassK0Short", "h1dMassK0Short", {HistType::kTH1F, {axisK0ShortMass}});
    histos.add("h1dMassLambda", "h1dMassLambda", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h1dMassAntiLambda", "h1dMassAntiLambda", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h1dMassCompatibleK0Short", "h1dMassCompatibleK0Short", {HistType::kTH1F, {axisK0ShortMass}});
    histos.add("h1dMassCompatibleLambda", "h1dMassCompatibleLambda", {HistType::kTH1F, {axisLambdaMass}});
    histos.add("h1dMassCompatibleAntiLambda", "h1dMassCompatibleAntiLambda", {HistType::kTH1F, {axisLambdaMass}});

    histos.add("h3dMassK0Short", "h3dMassK0Short", {HistType::kTH3F, {centAxis, axisPt, axisK0ShortMass}});
    histos.add("h3dMassLambda", "h3dMassLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});
    histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});

    histos.add("h3dMassCompatibleK0Short", "h3dMassCompatibleK0Short", {HistType::kTH3F, {centAxis, axisPt, axisK0ShortMass}});
    histos.add("h3dMassCompatibleLambda", "h3dMassCompatibleLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});
    histos.add("h3dMassCompatibleAntiLambda", "h3dMassCompatibleAntiLambda", {HistType::kTH3F, {centAxis, axisPt, axisLambdaMass}});

    // --- ASSOCIATED ---
    // V0 Radius
    if (doprocessSim) {
      histos.add("h2dAssocLambdaRadiusVsPt", "hLambdaRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});

      // Invariant Mass
      histos.add("h2dAssocLambdaMassVsPt", "hLambdaMassVsPt", {HistType::kTH2F, {axisPt, axisLambdaMass}});

      // radius vs prong length
      histos.add("h2dAssocProtonLengthVsRadius", "h2dAssocProtonLengthVsRadius", {HistType::kTH2F, {axisRadius, axisLength}});
      histos.add("h2dAssocPionLengthVsRadius", "h2dAssocPionLengthVsRadius", {HistType::kTH2F, {axisRadius, axisLength}});

      // recalculated vs topv lengths
      histos.add("h2dAssocProtonLengthVsLengthToPV", "h2dAssocProtonLengthVsRadiusToPV", {HistType::kTH2F, {axisRadius, axisRadius}});
      histos.add("h2dAssocPionLengthVsLengthToPV", "h2dAssocPionLengthVsLengthToPV", {HistType::kTH2F, {axisRadius, axisRadius}});

      // TOF PID testing for prongs
      histos.add("h2dAssocProtonDeltaTimeVsPt", "h2dAssocProtonDeltaTimeVsPt", {HistType::kTH2F, {axisPt, axisDeltaTime}});
      histos.add("h2dAssocProtonDeltaTimeVsRadius", "h2dAssocProtonDeltaTimeVsRadius", {HistType::kTH2F, {axisRadius, axisDeltaTime}});
      histos.add("h2dAssocPionDeltaTimeVsPt", "h2dAssocPionDeltaTimeVsPt", {HistType::kTH2F, {axisPt, axisDeltaTime}});
      histos.add("h2dAssocPionDeltaTimeVsRadius", "h2dAssocPionDeltaTimeVsRadius", {HistType::kTH2F, {axisRadius, axisDeltaTime}});

      // delta lambda decay time
      histos.add("h2dAssocLambdaDeltaDecayTime", "h2dAssocLambdaDeltaDecayTime", {HistType::kTH2F, {axisPt, axisDeltaTime}});
      histos.add("h2dAssocLambdaDeltaDecayTime_MassSelected", "h2dAssocLambdaDeltaDecayTime_MassSelected", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    }

    if (doprocessCascades) {
      histos.add("h1dMassXiMinus", "h1dMassXiMinus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassXiPlus", "h1dMassXiPlus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassOmegaMinus", "h1dMassOmegaMinus", {HistType::kTH1F, {massAxisOmega}});
      histos.add("h1dMassOmegaPlus", "h1dMassOmegaPlus", {HistType::kTH1F, {massAxisOmega}});
      histos.add("h1dMassCompatibleXiMinus", "h1dMassCompatibleXiMinus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassCompatibleXiPlus", "h1dMassCompatibleXiPlus", {HistType::kTH1F, {massAxisXi}});
      histos.add("h1dMassCompatibleOmegaMinus", "h1dMassCompatibleOmegaMinus", {HistType::kTH1F, {massAxisOmega}});
      histos.add("h1dMassCompatibleOmegaPlus", "h1dMassCompatibleOmegaPlus", {HistType::kTH1F, {massAxisOmega}});

      histos.add("h3dMassXiMinus", "h3dMassXiMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassXiPlus", "h3dMassXiPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassOmegaMinus", "h3dMassOmegaMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});
      histos.add("h3dMassOmegaPlus", "h3dMassOmegaPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});

      histos.add("h3dMassCompatibleXiMinus", "h3dMassCompatibleXiMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassCompatibleXiPlus", "h3dMassCompatibleXiPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisXi}});
      histos.add("h3dMassCompatibleOmegaMinus", "h3dMassCompatibleOmegaMinus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});
      histos.add("h3dMassCompatibleOmegaPlus", "h3dMassCompatibleOmegaPlus", {HistType::kTH3F, {centAxis, axisPt, massAxisOmega}});
    }
  }

  void processReal(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFBetas, aod::V0TOFDebugs, aod::V0TOFNSigmas> const& v0s, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {
    if (coll.centFT0C() > maxCentrality || coll.centFT0C() < minCentrality)
      return;

    for (auto& lambda : v0s) { // selecting photons from Sigma0

      if (TMath::Abs(lambda.eta()) > 0.5)
        continue;

      auto negExtra = lambda.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = lambda.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassLambda"), coll.centFT0C(), lambda.pt(), lambda.mLambda());
        histos.fill(HIST("h1dMassLambda"), lambda.mLambda());
        if (lambda.tofLambdaCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleLambda"), coll.centFT0C(), lambda.pt(), lambda.mLambda());
          histos.fill(HIST("h1dMassCompatibleLambda"), lambda.mLambda());
        }
      }

      if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassAntiLambda"), coll.centFT0C(), lambda.pt(), lambda.mAntiLambda());
        histos.fill(HIST("h1dMassAntiLambda"), lambda.mAntiLambda());
        if (lambda.tofAntiLambdaCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleAntiLambda"), coll.centFT0C(), lambda.pt(), lambda.mAntiLambda());
          histos.fill(HIST("h1dMassCompatibleAntiLambda"), lambda.mAntiLambda());
        }
      }

      if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaPion) {
        // lambda case
        histos.fill(HIST("h3dMassK0Short"), coll.centFT0C(), lambda.pt(), lambda.mK0Short());
        histos.fill(HIST("h1dMassK0Short"), lambda.mK0Short());
        if (lambda.tofK0ShortCompatibility(tofNsigmaCompatibility.value)) {
          histos.fill(HIST("h3dMassCompatibleK0Short"), coll.centFT0C(), lambda.pt(), lambda.mK0Short());
          histos.fill(HIST("h1dMassCompatibleK0Short"), lambda.mK0Short());
        }
      }

      histos.fill(HIST("h2dLambdaMassVsTOFCut"), lambda.mLambda(), TMath::Abs(lambda.posTOFDeltaTLaPr()));
      histos.fill(HIST("h2dLambdaMassVsTOFCutMeson"), lambda.mLambda(), TMath::Abs(lambda.negTOFDeltaTLaPi()));

      if (lambda.v0radius() < minV0Radius)
        continue;

      histos.fill(HIST("h2dLambdaDeltaDecayTime"), lambda.pt(), lambda.deltaDecayTimeLambda());
      if (TMath::Abs(lambda.mLambda() - 1.115683) < 0.01 && lambda.v0cosPA() > minCosPA) {
        histos.fill(HIST("h2dLambdaDeltaDecayTime_MassSelected"), lambda.pt(), lambda.deltaDecayTimeLambda());
      }

      if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
        histos.fill(HIST("hLambdaMass"), lambda.mLambda());

      histos.fill(HIST("h2dLambdaRadiusVsPt"), lambda.pt(), lambda.v0radius());
      histos.fill(HIST("h2dLambdaMassVsPt"), lambda.pt(), lambda.mLambda());

      histos.fill(HIST("h2dProtonDeltaTimeVsPt"), lambda.pt(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dProtonDeltaTimeVsRadius"), lambda.v0radius(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dPionDeltaTimeVsPt"), lambda.pt(), lambda.negTOFDeltaTLaPi());
      histos.fill(HIST("h2dPionDeltaTimeVsRadius"), lambda.v0radius(), lambda.negTOFDeltaTLaPi());
      histos.fill(HIST("h2dLambdaBeta"), lambda.p(), lambda.tofBetaLambda());

      if (TMath::Abs(lambda.mLambda() - 1.115683) < 0.01 && lambda.v0cosPA() > minCosPA) {
        histos.fill(HIST("h2dProtonDeltaTimeVsPt_MassSelected"), lambda.pt(), lambda.posTOFDeltaTLaPr());
        histos.fill(HIST("h2dProtonDeltaTimeVsRadius_MassSelected"), lambda.v0radius(), lambda.posTOFDeltaTLaPr());
        histos.fill(HIST("h2dPionDeltaTimeVsPt_MassSelected"), lambda.pt(), lambda.negTOFDeltaTLaPi());
        histos.fill(HIST("h2dPionDeltaTimeVsRadius_MassSelected"), lambda.v0radius(), lambda.negTOFDeltaTLaPi());
        histos.fill(HIST("h2dLambdaBeta_MassSelected"), lambda.p(), lambda.tofBetaLambda());
      }

      // Standard selection of time
      if (TMath::Abs(lambda.deltaDecayTimeLambda()) < maxDeltaTimeDecay) {
        histos.fill(HIST("h2dLambdaMassVsPt_DeltaDecayTime"), lambda.pt(), lambda.mLambda());
        if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
          histos.fill(HIST("hLambdaMass_DeltaDecayTime"), lambda.mLambda());
      }

      if (TMath::Abs(lambda.posTOFDeltaTLaPr()) < maxDeltaTimeProton) {
        histos.fill(HIST("h2dLambdaMassVsPt_ProtonTOF"), lambda.pt(), lambda.mLambda());
        if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
          histos.fill(HIST("hLambdaMass_ProtonTOF"), lambda.mLambda());
      }
      if (TMath::Abs(lambda.negTOFDeltaTLaPi()) < maxDeltaTimeProton) {
        histos.fill(HIST("h2dLambdaMassVsPt_PionTOF"), lambda.pt(), lambda.mLambda());
        if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
          histos.fill(HIST("hLambdaMass_PionTOF"), lambda.mLambda());
        if (TMath::Abs(lambda.posTOFDeltaTLaPr()) < maxDeltaTimeProton) {
          histos.fill(HIST("h2dLambdaMassVsPt_AllTOF"), lambda.pt(), lambda.mLambda());
          if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
            histos.fill(HIST("hLambdaMass_AllTOF"), lambda.mLambda());
        }
      }
      // Inversion of time selection for debug
      // if(TMath::Abs(lambda.deltaDecayTimeLambda())>maxDeltaTimeDecay && TMath::Abs(lambda.deltaDecayTimeLambda()) < 5e+4){
      //   histos.fill(HIST("h2dLambdaMassVsPt_DeltaDecayTime"), lambda.pt(), lambda.mLambda());
      //   if(lambda.pt()>minPtV0 && lambda.pt() < maxPtV0) histos.fill(HIST("hLambdaMass_DeltaDecayTime"), lambda.mLambda());
      // }

      if (TMath::Abs(lambda.posTOFDeltaTLaPr()) > maxDeltaTimeProton && TMath::Abs(lambda.posTOFDeltaTLaPr()) < 5e+4) {
        histos.fill(HIST("h2dLambdaMassVsPt_InvertProtonTOF"), lambda.pt(), lambda.mLambda());
        if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
          histos.fill(HIST("hLambdaMass_InvertProtonTOF"), lambda.mLambda());
      }
      if (TMath::Abs(lambda.negTOFDeltaTLaPi()) > maxDeltaTimeProton && TMath::Abs(lambda.negTOFDeltaTLaPi()) < 5e+4) {
        histos.fill(HIST("h2dLambdaMassVsPt_InvertPionTOF"), lambda.pt(), lambda.mLambda());
        if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
          histos.fill(HIST("hLambdaMass_InvertPionTOF"), lambda.mLambda());
        if (TMath::Abs(lambda.posTOFDeltaTLaPr()) > maxDeltaTimeProton && TMath::Abs(lambda.posTOFDeltaTLaPr()) < 5e+4) {
          histos.fill(HIST("h2dLambdaMassVsPt_InvertAllTOF"), lambda.pt(), lambda.mLambda());
          if (lambda.pt() > minPtV0 && lambda.pt() < maxPtV0)
            histos.fill(HIST("hLambdaMass_InvertAllTOF"), lambda.mLambda());
        }
      }
    }
  }

  void processSim(aod::StraCollision const&, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas, aod::V0TOFPIDs, aod::V0TOFBetas> const& v0s)
  {
    for (auto& lambda : v0s) { // selecting photons from Sigma0
      if (lambda.pdgCode() != 3122)
        continue;
      if (!lambda.isPhysicalPrimary())
        continue;
      if (lambda.pdgCodePositive() != 2212)
        continue;
      if (lambda.pdgCodeNegative() != -211)
        continue;

      histos.fill(HIST("hAssocLambdaMass"), lambda.mLambda());

      histos.fill(HIST("h2dAssocLambdaRadiusVsPt"), lambda.pt(), lambda.v0radius());
      histos.fill(HIST("h2dAssocLambdaMassVsPt"), lambda.pt(), lambda.mLambda());

      histos.fill(HIST("h2dAssocProtonDeltaTimeVsPt"), lambda.pt(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dAssocProtonDeltaTimeVsRadius"), lambda.v0radius(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dAssocPionDeltaTimeVsPt"), lambda.pt(), lambda.negTOFDeltaTLaPi());
      histos.fill(HIST("h2dAssocPionDeltaTimeVsRadius"), lambda.v0radius(), lambda.negTOFDeltaTLaPi());

      // delta lambda decay time
      histos.fill(HIST("h2dAssocLambdaDeltaDecayTime"), lambda.pt(), lambda.deltaDecayTimeLambda());
    }
  }

  // ____________________________________________________________________________
  // QA TOF NSigma quantities

  Partition<soa::Join<aod::CascCores, aod::CascExtras>> negCasc = aod::cascdata::sign < 0;
  Partition<soa::Join<aod::CascCores, aod::CascExtras>> posCasc = aod::cascdata::sign > 0;

  Filter preFilter =
    nabs(aod::cascdata::dcapostopv) > v0setting_dcapostopv&& nabs(aod::cascdata::dcanegtopv) > v0setting_dcanegtopv&& nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv&& aod::cascdata::dcaV0daughters < v0setting_dcav0dau&& aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau;

  void processCascades(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& col, soa::Filtered<soa::Join<aod::CascCores, aod::CascCollRefs, aod::CascExtras, aod::CascTOFNSigmas>> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {
    histos.fill(HIST("hEventCentrality"), col.centFT0C());
    if (col.centFT0C() > maxCentrality || col.centFT0C() < minCentrality)
      return;

    for (auto& casc : Cascades) {
      // major selections here
      if (casc.v0radius() > v0setting_radius &&
          casc.cascradius() > cascadesetting_cascradius &&
          casc.v0cosPA(col.posX(), col.posY(), col.posZ()) > v0setting_cospa &&
          casc.casccosPA(col.posX(), col.posY(), col.posZ()) > cascadesetting_cospa &&
          casc.dcav0topv(col.posX(), col.posY(), col.posZ()) > cascadesetting_mindcav0topv &&
          TMath::Abs(casc.mLambda() - 1.115683) < cascadesetting_v0masswindow) {

        auto negExtra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
        auto posExtra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
        auto bachExtra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

        if (negExtra.tpcCrossedRows() < tpcCrossedRows || posExtra.tpcCrossedRows() < tpcCrossedRows || bachExtra.tpcCrossedRows() < tpcCrossedRows)
          continue;

        if (casc.sign() < 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(bachExtra.tpcNSigmaPi()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassXiMinus"), col.centFT0C(), casc.pt(), casc.mXi());
            histos.fill(HIST("h1dMassXiMinus"), casc.mXi());
            if (casc.tofXiCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleXiMinus"), col.centFT0C(), casc.pt(), casc.mXi());
              histos.fill(HIST("h1dMassCompatibleXiMinus"), casc.mXi());
            }
          }
          if (TMath::Abs(posExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(negExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(bachExtra.tpcNSigmaKa()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassOmegaMinus"), col.centFT0C(), casc.pt(), casc.mOmega());
            histos.fill(HIST("h1dMassOmegaMinus"), casc.mOmega());
            if (casc.tofOmegaCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleOmegaMinus"), col.centFT0C(), casc.pt(), casc.mOmega());
              histos.fill(HIST("h1dMassCompatibleOmegaMinus"), casc.mOmega());
            }
          }
        } else {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(bachExtra.tpcNSigmaPi()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassXiPlus"), col.centFT0C(), casc.pt(), casc.mXi());
            histos.fill(HIST("h1dMassXiPlus"), casc.mXi());
            if (casc.tofXiCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleXiPlus"), col.centFT0C(), casc.pt(), casc.mXi());
              histos.fill(HIST("h1dMassCompatibleXiPlus"), casc.mXi());
            }
          }

          if (TMath::Abs(posExtra.tpcNSigmaPi()) < tpcNsigmaPion && TMath::Abs(negExtra.tpcNSigmaPr()) < tpcNsigmaProton && TMath::Abs(bachExtra.tpcNSigmaKa()) < tpcNsigmaBachelor) {
            histos.fill(HIST("h3dMassOmegaPlus"), col.centFT0C(), casc.pt(), casc.mOmega());
            histos.fill(HIST("h1dMassOmegaPlus"), casc.mOmega());
            if (casc.tofOmegaCompatibility(tofNsigmaCompatibilityCascades.value)) {
              histos.fill(HIST("h3dMassCompatibleOmegaPlus"), col.centFT0C(), casc.pt(), casc.mOmega());
              histos.fill(HIST("h1dMassCompatibleOmegaPlus"), casc.mOmega());
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(strangepidqa, processReal, "Produce real information", true);
  PROCESS_SWITCH(strangepidqa, processSim, "Produce simulated information", true);
  PROCESS_SWITCH(strangepidqa, processCascades, "Process real cascades", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<strangepidqa>(cfgc)};
}
