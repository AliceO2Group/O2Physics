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

  // Invariant Mass
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

  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});

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
  }

  void processReal(soa::Join<aod::StraCollisions, aod::StraCents>::iterator const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFs, aod::V0TOFPIDs, aod::V0TOFBetas, aod::V0TOFDebugs> const& v0s)
  {
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    if (coll.centFT0C() > maxCentrality || coll.centFT0C() < minCentrality)
      return;

    for (auto& lambda : v0s) { // selecting photons from Sigma0

      if (TMath::Abs(lambda.eta()) > 0.5)
        continue;

      if (TMath::Abs(lambda.mLambda() - 1.115683) < 0.01 && lambda.v0cosPA() > minCosPA) {
        histos.fill(HIST("hPositiveStatus"), 0.0f);
        if (lambda.posTOFSignal() > 1e-5 && lambda.posTOFEventTime() < 0.0f)
          histos.fill(HIST("hPositiveStatus"), 1.0f);
        if (lambda.posTOFSignal() < 0.0f && lambda.posTOFEventTime() > 1e-5)
          histos.fill(HIST("hPositiveStatus"), 2.0f);
        if (lambda.posTOFSignal() > 1e-5 && lambda.posTOFEventTime() > 1e-5)
          histos.fill(HIST("hPositiveStatus"), 3.0f);

        histos.fill(HIST("hNegativeStatus"), 0.0f);
        if (lambda.negTOFSignal() > 1e-5 && lambda.negTOFEventTime() < 0.0f)
          histos.fill(HIST("hNegativeStatus"), 1.0f);
        if (lambda.negTOFSignal() < 0.0f && lambda.negTOFEventTime() > 1e-5)
          histos.fill(HIST("hNegativeStatus"), 2.0f);
        if (lambda.negTOFSignal() > 1e-5 && lambda.negTOFEventTime() > 1e-5)
          histos.fill(HIST("hNegativeStatus"), 3.0f);

        if (lambda.posTOFLength() > 0.0f) {
          histos.fill(HIST("hPositiveStatusReachedTOF"), 0.0f);
          if (lambda.posTOFSignal() > 1e-5 && lambda.posTOFEventTime() < 0.0f)
            histos.fill(HIST("hPositiveStatusReachedTOF"), 1.0f);
          if (lambda.posTOFSignal() < 0.0f && lambda.posTOFEventTime() > 1e-5)
            histos.fill(HIST("hPositiveStatusReachedTOF"), 2.0f);
          if (lambda.posTOFSignal() > 1e-5 && lambda.posTOFEventTime() > 1e-5)
            histos.fill(HIST("hPositiveStatusReachedTOF"), 3.0f);
        }

        if (lambda.negTOFLength() > 0.0f) {
          histos.fill(HIST("hNegativeStatusReachedTOF"), 0.0f);
          if (lambda.negTOFSignal() > 1e-5 && lambda.negTOFEventTime() < 0.0f)
            histos.fill(HIST("hNegativeStatusReachedTOF"), 1.0f);
          if (lambda.negTOFSignal() < 0.0f && lambda.negTOFEventTime() > 1e-5)
            histos.fill(HIST("hNegativeStatusReachedTOF"), 2.0f);
          if (lambda.negTOFSignal() > 1e-5 && lambda.negTOFEventTime() > 1e-5)
            histos.fill(HIST("hNegativeStatusReachedTOF"), 3.0f);
        }

        // properties of the positive and negative prongs / debug
        histos.fill(HIST("h2dTimeVsLengthProtonProng"), lambda.posTOFLength(), lambda.posLifetimePr());
        histos.fill(HIST("h2dTimeVsLengthPionProng"), lambda.negTOFLength(), lambda.negLifetimePr());
      }

      histos.fill(HIST("h2dLambdaMassVsTOFCut"), lambda.mLambda(), TMath::Abs(lambda.posTOFDeltaTLaPr()));
      histos.fill(HIST("h2dLambdaMassVsTOFCutMeson"), lambda.mLambda(), TMath::Abs(lambda.negTOFDeltaTLaPi()));

      if (lambda.v0radius() < minV0Radius)
        continue;
      if (lambda.posTOFSignal() < 0 && requireTOFsignalProton)
        continue;
      if (lambda.negTOFSignal() < 0 && requireTOFsignalPion)
        continue;

      histos.fill(HIST("h2dLambdaDeltaDecayTime"), lambda.pt(), lambda.deltaDecayTimeLambda());
      if (TMath::Abs(lambda.mLambda() - 1.115683) < 0.01 && lambda.v0cosPA() > minCosPA) {
        histos.fill(HIST("h2dLambdaDeltaDecayTime_MassSelected"), lambda.pt(), lambda.deltaDecayTimeLambda());
      }

      if (lambda.posTOFEventTime() < 0 && requireTOFEventTimeProton)
        continue;
      if (lambda.negTOFEventTime() < 0 && requireTOFEventTimePion)
        continue;

      histos.fill(HIST("h2dLambdaMassVsTOFCutWithSignals"), lambda.mLambda(), TMath::Abs(lambda.posTOFDeltaTLaPr()));
      histos.fill(HIST("h2dLambdaMassVsTOFCutMesonWithSignals"), lambda.mLambda(), TMath::Abs(lambda.negTOFDeltaTLaPi()));

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

      histos.fill(HIST("h2dProtonLengthVsRadius"), lambda.v0radius(), lambda.posTOFLength());
      histos.fill(HIST("h2dPionLengthVsRadius"), lambda.v0radius(), lambda.negTOFLength());

      histos.fill(HIST("h2dProtonLengthVsLengthToPV"), lambda.posTOFLengthToPV(), lambda.posTOFLength());
      histos.fill(HIST("h2dPionLengthVsLengthToPV"), lambda.negTOFLengthToPV(), lambda.negTOFLength());

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

  void processSim(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas, aod::V0TOFs, aod::V0TOFPIDs, aod::V0TOFBetas> const& v0s)
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

      if (TMath::Abs(lambda.negTOFDeltaTLaPi()) > 800) {
        histos.fill(HIST("hAssocLambdaMassGoodTime"), lambda.mLambda());
      } else {
        histos.fill(HIST("hAssocLambdaMassBadTime"), lambda.mLambda());
      }

      histos.fill(HIST("h2dAssocLambdaRadiusVsPt"), lambda.pt(), lambda.v0radius());
      histos.fill(HIST("h2dAssocLambdaMassVsPt"), lambda.pt(), lambda.mLambda());

      histos.fill(HIST("h2dAssocProtonDeltaTimeVsPt"), lambda.pt(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dAssocProtonDeltaTimeVsRadius"), lambda.v0radius(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dAssocPionDeltaTimeVsPt"), lambda.pt(), lambda.negTOFDeltaTLaPi());
      histos.fill(HIST("h2dAssocPionDeltaTimeVsRadius"), lambda.v0radius(), lambda.negTOFDeltaTLaPi());

      histos.fill(HIST("h2dAssocProtonLengthVsRadius"), lambda.v0radius(), lambda.posTOFLength());
      histos.fill(HIST("h2dAssocPionLengthVsRadius"), lambda.v0radius(), lambda.negTOFLength());

      histos.fill(HIST("h2dAssocProtonLengthVsLengthToPV"), lambda.posTOFLengthToPV(), lambda.posTOFLength());
      histos.fill(HIST("h2dAssocPionLengthVsLengthToPV"), lambda.negTOFLengthToPV(), lambda.negTOFLength());

      // delta lambda decay time
      histos.fill(HIST("h2dAssocLambdaDeltaDecayTime"), lambda.pt(), lambda.deltaDecayTimeLambda());
    }
  }

  PROCESS_SWITCH(strangepidqa, processReal, "Produce real information", true);
  PROCESS_SWITCH(strangepidqa, processSim, "Produce simulated information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<strangepidqa>(cfgc)};
}
