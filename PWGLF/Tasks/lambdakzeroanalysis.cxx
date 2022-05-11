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
// Example V0 analysis task
// ========================
//
// This code loops over a V0Data table and produces some
// standard analysis output. It requires either
// the lambdakzerofinder or the lambdakzeroproducer tasks
// to have been executed in the workflow (before).
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "PID/PIDResponse.h"

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

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPr>;

struct lambdakzeroQa {
  //Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hArmenterosPreAnalyserCuts", "hArmenterosPreAnalyserCuts", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}}},
    },
  };
  void init(InitContext const&)
  {
    AxisSpec massAxisK0Short = {600, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"};
    AxisSpec massAxisLambda = {600, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {massAxisK0Short}});
    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {massAxisLambda}});
    registry.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {massAxisLambda}});
  }
  void process(aod::Collision const& collision, aod::V0Datas const& fullV0s)
  {

    for (auto& v0 : fullV0s) {
      registry.fill(HIST("hMassK0Short"), v0.mK0Short());
      registry.fill(HIST("hMassLambda"), v0.mLambda());
      registry.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());

      registry.fill(HIST("hV0Radius"), v0.v0radius());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
      registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());
      registry.fill(HIST("hArmenterosPreAnalyserCuts"), v0.alpha(), v0.qtarm());
    }
  }
};

struct lambdakzeroAnalysis {

  HistogramRegistry registry{
    "registry",
    {
      {"h3dMassK0Short", "h3dMassK0Short", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 0.450f, 0.550f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dMassLambda", "h3dMassLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 1.015f, 1.215f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dMassAntiLambda", "h3dMassAntiLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 1.015f, 1.215f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hSelectedEventCounter", "hSelectedEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      {"hArmenterosPostAnalyserCuts", "hArmenterosPostAnalyserCuts", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}}},
    },
  };

  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};

  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisK0Short = {200, 0.450f, 0.550f, "Inv. Mass (GeV/c^{2})"};
    AxisSpec massAxisLambda = {200, 1.015f, 1.215f, "Inv. Mass (GeV/c^{2})"};

    registry.add("h3dMassK0ShortDca", "h3dMassK0ShortDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisK0Short}});
    registry.add("h3dMassLambdaDca", "h3dMassLambdaDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});
    registry.add("h3dMassAntiLambdaDca", "h3dMassAntiLambdaDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});
  }

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<int> saveDcaHist{"saveDcaHist", 0, "saveDcaHist"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<bool> boolArmenterosCut{"boolArmenterosCut", true, "cut on Armenteros-Podolanski graph"};
  Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  // void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s) //for now CentV0M info is not available for run 3 pp
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, MyTracks const& tracks)
  {
    if (eventSelection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hSelectedEventCounter"), 0.5);

    for (auto& v0 : fullV0s) {
      //FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        if (TMath::Abs(v0.yLambda()) < rapidity) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kLambda0) < lifetimecut->get("lifetimecutLambda")) {

            // Lambda
            if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPr()) < TpcPidNsigmaCut) { //previous 900Gev pp analysis had nSigma< 5 for pt<0.7Gev and tpcNSigmaStorePi<3 for pt>0.7GeV; and no cut on K0S
              registry.fill(HIST("h3dMassLambda"), 0., v0.pt(), v0.mLambda());            //collision.centV0M() instead of 0. once available
              registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mLambda());
              }
            }

            // AntiLambda
            if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPr()) < TpcPidNsigmaCut) { //previous 900Gev pp analysis had nSigma< 5 for pt<0.7Gev and tpcNSigmaStorePi<3 for pt>0.7GeV; and no cut on K0S
              registry.fill(HIST("h3dMassAntiLambda"), 0., v0.pt(), v0.mAntiLambda());
              registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassAntiLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiLambda());
              }
            }
          }
        }

        //K0Short
        if (TMath::Abs(v0.yK0Short()) < rapidity) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kK0Short) < lifetimecut->get("lifetimecutK0S")) {
            if ((v0.qtarm() > paramArmenterosCut * TMath::Abs(v0.alpha())) || !boolArmenterosCut) {
              registry.fill(HIST("h3dMassK0Short"), 0., v0.pt(), v0.mK0Short());
              registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassK0ShortDca"), v0.dcaV0daughters(), v0.pt(), v0.mK0Short());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(lambdakzeroAnalysis, processRun3, "Process Run 3 data", true);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, MyTracks const& tracks)
  {
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    registry.fill(HIST("hSelectedEventCounter"), 0.5);

    for (auto& v0 : fullV0s) {
      //FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        if (TMath::Abs(v0.yLambda()) < rapidity) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kLambda0) < lifetimecut->get("lifetimecutLambda")) {

            // Lambda
            if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPr()) < TpcPidNsigmaCut) { //previous 900Gev pp analysis had nSigma< 5 for pt<0.7Gev and tpcNSigmaStorePi<3 for pt>0.7GeV; and no cut on K0S
              registry.fill(HIST("h3dMassLambda"), collision.centRun2V0M(), v0.pt(), v0.mLambda());
              registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mLambda());
              }
            }

            // AntiLambda
            if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPr()) < TpcPidNsigmaCut) { //previous 900Gev pp analysis had nSigma< 5 for pt<0.7Gev and tpcNSigmaStorePi<3 for pt>0.7GeV; and no cut on K0S
              registry.fill(HIST("h3dMassAntiLambda"), collision.centRun2V0M(), v0.pt(), v0.mAntiLambda());
              registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassAntiLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiLambda());
              }
            }
          }
        }

        //K0Short
        if (TMath::Abs(v0.yK0Short()) < rapidity) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kK0Short) < lifetimecut->get("lifetimecutK0S")) {
            if ((v0.qtarm() > paramArmenterosCut * v0.alpha()) || !boolArmenterosCut) {
              registry.fill(HIST("h3dMassK0Short"), collision.centRun2V0M(), v0.pt(), v0.mK0Short());
              registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassK0ShortDca"), v0.dcaV0daughters(), v0.pt(), v0.mK0Short());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(lambdakzeroAnalysis, processRun2, "Process Run 2 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeroAnalysis>(cfgc),
    adaptAnalysisTask<lambdakzeroQa>(cfgc)};
}
