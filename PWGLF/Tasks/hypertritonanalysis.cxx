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
// the hypertritonfinder or the hypertritonproducer tasks
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

//using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCHe, aod::pidTPCTr, aod::pidTPCKa, aod::pidTPCPr>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTPCFullKa, aod::pidTPCFullPr>;

struct hypertritonQa {
  //Basic checks
  HistogramRegistry registry{
    "registry",
      {
        {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
        {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
        {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
        {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, 10.0f, 10.0f, "cm"}}}},
        {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
        {"hArmenterosPreAnalyserCuts", "hArmenterosPreAnalyserCuts", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}}},
      },
  };
  void init(InitContext const&)
  {
    AxisSpec massAxis = {300, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hMassCandidates", "hMassCandidates", {HistType::kTH1F, {massAxis}});
  }
  void process(aod::Collision const& collision, aod::V0Datas const& fullV0s)
  {

    for (auto& v0 : fullV0s) {
      registry.fill(HIST("hV0Radius"), v0.v0radius());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
      registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());
      registry.fill(HIST("hArmenterosPreAnalyserCuts"), v0.alpha(), v0.qtarm());
      registry.fill(HIST("hMassCandidates"), v0.mHypertriton());
    }
  }
};

struct hypertritonAnalysis {

  HistogramRegistry registry{
    "registry",
      {
        {"hSelectedEventCounter", "hSelectedEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
        {"hSelectedV0Counter", "hSelectedV0Counter", {HistType::kTH1F, {{8, 0.0f, 8.0f}}}},
        {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
        {"hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
        {"hMassHypertritonTotal", "hMassHypertritonTotal", {HistType::kTH1F, {{160, 2.8f, 3.2f}}}},
        {"hNSigmaHelium3", "hNSigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaPion", "hNSigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaTriton", "hNSigmaTriton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaKaon", "hNSigmaKaon", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaProton", "hNSigmaProton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hPtHelium3", "hPtHelium3", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hPtPion", "hPtPion", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hPtAntiHelium3", "hPtAntiHelium3", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hPtAntiPion", "hPtAntiPion", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"h3dMassHypertriton", "h3dMassHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"h3dMassAntiHypertriton", "h3dMassAntiHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"hArmenterosPostAnalyserCuts", "hArmenterosPostAnalyserCuts", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}}},
      },
  };

  Configurable<int> saveDcaHist{"saveDcaHist", 0, "saveDcaHist"};

  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};

  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisHypertriton = {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"};

    if (saveDcaHist==1){
      registry.add("h3dMassHypertritonDca", "h3dMassHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
      registry.add("h3dMassAntiHypertritonDca", "h3dMassAntiHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
    }

    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(2, "V0Radius");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(3, "V0CosPA");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(4, "Rapidity");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(5, "Lifetime");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(6, "DcaV0Dau");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(7, "TPCPID");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(8, "PionDcatoPV");
  }

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};//loose cut
  Configurable<float> dcapiontopv{"dcapiontopv", .1, "DCA Pion To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.9, "rapidity"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  //Configurable<bool> boolArmenterosCut{"boolArmenterosCut", true, "cut on Armenteros-Podolanski graph"};//unknown
  //Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};//unknown
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<float> lifetimecut{"lifetimecut", 25., "lifetimecut"}; //ct

  //Filter dcaFilterV0 = aod::v0data::dcaV0daughters < dcav0dau;

  //Dynamic columns doesn't work in filters
  //Filter cpaFilterV0 = aod::v0data::v0cosPA > v0cospa;
  //Filter radiusFilterV0 = aod::v0data::v0radius > v0radius;
  //Filter YFilterV0 = nabs(aod::v0data::yHypertriton) < rapidity;// need to be checked
  //Pt Cut; TPCTritonNSigmaCut; Pion PID;

  // void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s) //for now CentV0M info is not available for run 3 pp
  //void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, MyTracks const& tracks)
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& fullV0s, MyTracks const& tracks)
  {
    registry.fill(HIST("hSelectedEventCounter"), 0.5);
    /*if (eventSelection && !collision.sel8()) {
      return;
      }
      registry.fill(HIST("hSelectedEventCounter"), 1.5);*/

    for (auto& v0 : fullV0s) {
      //FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      registry.fill(HIST("hSelectedV0Counter"), 0.5);
      if (v0.v0radius() < v0radius){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 1.5);
      if(v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 2.5);
      if (TMath::Abs(v0.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 3.5);
      if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991 > lifetimecut) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 4.5);
      if (v0.dcaV0daughters() > dcav0dau){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 5.5);

      registry.fill(HIST("hNSigmaHelium3"), v0.posTrack_as<MyTracks>().tpcNSigmaHe());
      registry.fill(HIST("hNSigmaHelium3"), v0.negTrack_as<MyTracks>().tpcNSigmaHe());
      registry.fill(HIST("hNSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaTriton"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaTriton"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaKaon"), v0.posTrack_as<MyTracks>().tpcNSigmaKa());
      registry.fill(HIST("hNSigmaKaon"), v0.negTrack_as<MyTracks>().tpcNSigmaKa());
      registry.fill(HIST("hNSigmaProton"), v0.posTrack_as<MyTracks>().tpcNSigmaPr());
      registry.fill(HIST("hNSigmaProton"), v0.negTrack_as<MyTracks>().tpcNSigmaPr());
      // Hypertriton
      if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
        registry.fill(HIST("hSelectedV0Counter"), 6.5);

        if (TMath::Abs(v0.dcanegtopv()) > dcapiontopv) {
          registry.fill(HIST("hSelectedV0Counter"), 7.5);

          registry.fill(HIST("hPtHelium3"), v0.positivept());
          registry.fill(HIST("hPtAntiPion"), v0.negativept());
          registry.fill(HIST("hMassHypertriton"), v0.mHypertriton());
          registry.fill(HIST("hMassHypertritonTotal"), v0.mHypertriton());
          registry.fill(HIST("h3dMassHypertriton"), 0., v0.pt(), v0.mHypertriton());            //collision.centV0M() instead of 0. once available
          registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
          if (saveDcaHist == 1) {
            registry.fill(HIST("h3dMassHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mHypertriton());
          }
        }
      }

      // AntiHypertriton
      if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {

        registry.fill(HIST("hSelectedV0Counter"), 6.5);
        if (TMath::Abs(v0.dcapostopv()) > dcapiontopv) {
          registry.fill(HIST("hSelectedV0Counter"), 7.5);

          registry.fill(HIST("hPtAntiHelium3"), v0.negativept());
          registry.fill(HIST("hPtPion"), v0.positivept());
          registry.fill(HIST("hMassAntiHypertriton"), v0.mAntiHypertriton());
          registry.fill(HIST("hMassHypertritonTotal"), v0.mAntiHypertriton());
          registry.fill(HIST("h3dMassAntiHypertriton"), 0., v0.pt(), v0.mAntiHypertriton());
          registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
          if (saveDcaHist == 1) {
            registry.fill(HIST("h3dMassAntiHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiHypertriton());
          }
        }
      }

    }
  }

  PROCESS_SWITCH(hypertritonAnalysis, processRun3, "Process Run 3 data", true);

  //void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, MyTracks const& tracks)
  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, aod::V0Datas const& fullV0s, MyTracks const& tracks)
  {
    registry.fill(HIST("hSelectedEventCounter"), 0.5);
    /*if (!collision.alias()[kINT7]) {
      return;
      }
      if (eventSelection && !collision.sel7()) {
      return;
      }
      registry.fill(HIST("hSelectedEventCounter"), 1.5);*/

    for (auto& v0 : fullV0s) {
      //FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      registry.fill(HIST("hSelectedV0Counter"), 0.5);
      if (v0.v0radius() < v0radius){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 1.5);
      if(v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 2.5);
      if (TMath::Abs(v0.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 3.5);
      if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991 > lifetimecut) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 4.5);
      if (v0.dcaV0daughters() > dcav0dau){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 5.5);

      registry.fill(HIST("hNSigmaHelium3"), v0.posTrack_as<MyTracks>().tpcNSigmaHe());
      registry.fill(HIST("hNSigmaHelium3"), v0.negTrack_as<MyTracks>().tpcNSigmaHe());
      registry.fill(HIST("hNSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaTriton"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaTriton"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());
      // Hypertriton
      if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
        registry.fill(HIST("hSelectedV0Counter"), 6.5);

        if (TMath::Abs(v0.dcanegtopv()) > dcapiontopv) {
          registry.fill(HIST("hSelectedV0Counter"), 7.5);

          registry.fill(HIST("hPtHelium3"), v0.positivept());
          registry.fill(HIST("hPtAntiPion"), v0.negativept());
          registry.fill(HIST("hMassHypertriton"), v0.mHypertriton());
          registry.fill(HIST("hMassHypertritonTotal"), v0.mHypertriton());
          registry.fill(HIST("h3dMassHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mHypertriton());
          registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
          if (saveDcaHist == 1) {
            registry.fill(HIST("h3dMassHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mHypertriton());
          }
        }
      }

      // AntiHypertriton
      if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {

        registry.fill(HIST("hSelectedV0Counter"), 6.5);
        if (TMath::Abs(v0.dcapostopv()) > dcapiontopv) {
          registry.fill(HIST("hSelectedV0Counter"), 7.5);

          registry.fill(HIST("hPtAntiHelium3"), v0.negativept());
          registry.fill(HIST("hPtPion"), v0.positivept());
          registry.fill(HIST("hMassAntiHypertriton"), v0.mAntiHypertriton());
          registry.fill(HIST("hMassHypertritonTotal"), v0.mAntiHypertriton());
          registry.fill(HIST("h3dMassAntiHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mAntiHypertriton());
          registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
          if (saveDcaHist == 1) {
            registry.fill(HIST("h3dMassAntiHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiHypertriton());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(hypertritonAnalysis, processRun2, "Process Run 2 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertritonAnalysis>(cfgc),
      adaptAnalysisTask<hypertritonQa>(cfgc)};
}
