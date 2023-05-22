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
// the hypertritonbuilder or hypertritonfinder tasks
// to have been executed in the workflow (before).
//
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include "DataFormatsTPC/BetheBlochAleph.h"

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

//using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCHe, aod::pidTPCTr, aod::pidTPCKa, aod::pidTPCPr>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTPCFullKa, aod::pidTPCFullPr>;

inline float GetTPCNSigmaHe3(float p, float TPCSignal)
{
  float bg = p/2.80839;
  return  (TPCSignal - o2::tpc::BetheBlochAleph(bg, -9.973f, -18.5543f, 29.5704f, 2.02064f, -3.85076f)) / (TPCSignal*0.0812);
}

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
        {"hV0Pt", "hV0Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
        {"hPosPt", "hPosPt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
        {"hNegPt", "hNegPt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
        {"hHelium3Pt", "hHelium3Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      },
  };
  void init(InitContext const&)
  {
    AxisSpec massAxis = {120, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {massAxis}});
    registry.add("hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {massAxis}});
  }
  void process(aod::Collision const& collision, aod::V0Datas const& fullV0s, MyTracks const& tracks)
  {

    for (auto& v0 : fullV0s) {
      registry.fill(HIST("hV0Radius"), v0.v0radius());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
      registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());
      registry.fill(HIST("hArmenterosPreAnalyserCuts"), v0.alpha(), v0.qtarm());
      registry.fill(HIST("hV0Pt"), v0.pt());
      registry.fill(HIST("hPosPt"), v0.positivept());
      registry.fill(HIST("hNegPt"), v0.negativept());
      if (TMath::Abs( GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal()) ) < 5){
        registry.fill(HIST("hMassHypertriton"), v0.mHypertriton());
        registry.fill(HIST("hHelium3Pt"), v0.positivept());
      }
      if (TMath::Abs( GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal()) ) < 5){
        registry.fill(HIST("hMassAntiHypertriton"), v0.mAntiHypertriton());
        registry.fill(HIST("hHelium3Pt"), v0.negativept());
      }
    }
  }
};

struct hypertritonAnalysis {

  HistogramRegistry registry{
    "registry",
      {
        {"hSelectedEventCounter", "hSelectedEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
        {"hSelectedV0Counter", "hSelectedV0Counter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
        {"hTestCounter", "hTestCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
        //{"hSelectedCandidatesCounter", "hSelectedCandidatesCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},//delete cut which are applied in v0 Selection
        {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
        {"hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
        {"hMassHypertritonTotal", "hMassHypertritonTotal", {HistType::kTH1F, {{120, 2.9f, 3.2f}}}},
        {"hNSigmaHelium3", "hNSigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaPion", "hNSigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaTriton", "hNSigmaTriton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaKaon", "hNSigmaKaon", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaProton", "hNSigmaProton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hPtHelium3", "hPtHelium3", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hPtPion", "hPtPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hPtAntiHelium3", "hPtAntiHelium3", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hPtAntiPion", "hPtAntiPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"h3dMassHypertriton", "h3dMassHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"h3dMassAntiHypertriton", "h3dMassAntiHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"h3dTotalHypertriton", "h3dTotalHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
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
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(2, "V0CosPA");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(3, "TrackEta");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(4, "MomRapidity");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(5, "Lifetime");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(6, "DcaV0Dau");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(7, "TPCPID");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(8, "PtCut");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(9, "PionDcatoPV");

    /*registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(1, "Readin");
      registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(2, "Rapidity");
      registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(3, "Lifetime");
      registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(4, "TPC PID");
      registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(5, "PionPt");
      registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(6, "PionDcatoPv");*/
  }

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};//loose cut
  Configurable<float> dcapiontopv{"dcapiontopv", .1, "DCA Pion To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> etacut{"etacut", 0.9, "etacut"};
  Configurable<float> rapidity{"rapidity", 0.8, "rapidity"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  //Configurable<bool> boolArmenterosCut{"boolArmenterosCut", true, "cut on Armenteros-Podolanski graph"};//unknown
  //Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};//unknown
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<float> lifetimecut{"lifetimecut", 40., "lifetimecut"}; //ct

  //Filter dcaFilterV0 = aod::v0data::dcaV0daughters < dcav0dau;

  //Dynamic columns doesn't work in filters
  //Filter cpaFilterV0 = aod::v0data::v0cosPA > v0cospa;
  //Filter radiusFilterV0 = aod::v0data::v0radius > v0radius;
  //Filter RapidityFilterV0 = nabs(aod::v0data::yHypertriton) < rapidity;// need to be checked

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
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 1.5);
      if ( TMath::Abs(v0.posTrack_as<MyTracks>().eta()) > etacut || TMath::Abs(v0.negTrack_as<MyTracks>().eta()) > etacut){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 2.5);
      if (TMath::Abs(v0.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 3.5);
      double ct = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991; 
      if (ct > lifetimecut) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 4.5);
      if (v0.dcaV0daughters() > dcav0dau){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 5.5);

      registry.fill(HIST("hNSigmaHelium3"), GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ) );
      registry.fill(HIST("hNSigmaHelium3"), GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) );
      registry.fill(HIST("hNSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaTriton"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaTriton"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaKaon"), v0.posTrack_as<MyTracks>().tpcNSigmaKa());
      registry.fill(HIST("hNSigmaKaon"), v0.negTrack_as<MyTracks>().tpcNSigmaKa());
      registry.fill(HIST("hNSigmaProton"), v0.posTrack_as<MyTracks>().tpcNSigmaPr());
      registry.fill(HIST("hNSigmaProton"), v0.negTrack_as<MyTracks>().tpcNSigmaPr());
      // Hypertriton
      if (TMath::Abs( GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ) ) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
        registry.fill(HIST("hSelectedV0Counter"), 6.5);

        registry.fill(HIST("hTestCounter"), 0.5);
        if(v0.negativept() > 0.2 && v0.negativept() < 1.2 ){
          registry.fill(HIST("hTestCounter"), 1.5);
          if (v0.positivept() > 1.8 && v0.positivept() < 10){
            registry.fill(HIST("hTestCounter"), 2.5);
            if (v0.pt() > 2 && v0.pt() < 9 ){
              registry.fill(HIST("hTestCounter"), 3.5);
            }
          }
        }

        if(v0.negativept() > 0.2 && v0.negativept() < 1.2 && v0.positivept() > 1.8 && v0.positivept() < 10 && v0.pt() > 2 && v0.pt() < 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 7.5);

          if (TMath::Abs(v0.dcanegtopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 8.5);

            registry.fill(HIST("hPtHelium3"), v0.positivept());
            registry.fill(HIST("hPtAntiPion"), v0.negativept());
            registry.fill(HIST("hMassHypertriton"), v0.mHypertriton());
            registry.fill(HIST("hMassHypertritonTotal"), v0.mHypertriton());
            registry.fill(HIST("h3dMassHypertriton"), 0., v0.pt(), v0.mHypertriton());            //collision.centV0M() instead of 0. once available
            registry.fill(HIST("h3dTotalHypertriton"), ct, v0.pt(), v0.mHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mHypertriton());
            }
          }
        }
      }

      // AntiHypertriton
      if (TMath::Abs( GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) ) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {

        registry.fill(HIST("hSelectedV0Counter"), 6.5);

        registry.fill(HIST("hTestCounter"), 0.5);
        if(v0.positivept() > 0.2 && v0.positivept() < 1.2 ){
          registry.fill(HIST("hTestCounter"), 1.5);
          if (v0.negativept() > 1.8 && v0.negativept() < 10){
            registry.fill(HIST("hTestCounter"), 2.5);
            if (v0.pt() > 2 && v0.pt() < 9 ){
              registry.fill(HIST("hTestCounter"), 3.5);
            }
          }
        }

        if(v0.positivept() > 0.2 && v0.positivept() < 1.2 && v0.negativept() > 1.8 && v0.negativept() < 10 && v0.pt() > 2 && v0.pt() < 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 7.5);
          if (TMath::Abs(v0.dcapostopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 8.5);

            registry.fill(HIST("hPtAntiHelium3"), v0.negativept());
            registry.fill(HIST("hPtPion"), v0.positivept());
            registry.fill(HIST("hMassAntiHypertriton"), v0.mAntiHypertriton());
            registry.fill(HIST("hMassHypertritonTotal"), v0.mAntiHypertriton());
            registry.fill(HIST("h3dMassAntiHypertriton"), 0., v0.pt(), v0.mAntiHypertriton());
            registry.fill(HIST("h3dTotalHypertriton"), ct, v0.pt(), v0.mAntiHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassAntiHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiHypertriton());
            }
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
      if(v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 1.5);
      if ( TMath::Abs(v0.posTrack_as<MyTracks>().eta()) > etacut || TMath::Abs(v0.negTrack_as<MyTracks>().eta()) > etacut){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 2.5);
      if (TMath::Abs(v0.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 3.5);
      double ct = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991; 
      if (ct > lifetimecut) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 4.5);
      if (v0.dcaV0daughters() > dcav0dau){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 5.5);

      registry.fill(HIST("hNSigmaHelium3"), GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ) );
      registry.fill(HIST("hNSigmaHelium3"), GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) );
      registry.fill(HIST("hNSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaTriton"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaTriton"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());
      // Hypertriton
      if (TMath::Abs( GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ) ) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
        registry.fill(HIST("hSelectedV0Counter"), 6.5);

        if(v0.negTrack_as<MyTracks>().pt() > 0.2 && v0.negTrack_as<MyTracks>().pt() < 1.2 && v0.posTrack_as<MyTracks>().pt() > 1.8 && v0.posTrack_as<MyTracks>().pt() < 10 && v0.pt() > 2 && v0.pt() < 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 7.5);
          if (TMath::Abs(v0.dcanegtopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 8.5);

            registry.fill(HIST("hPtHelium3"), v0.positivept());
            registry.fill(HIST("hPtAntiPion"), v0.negativept());
            registry.fill(HIST("hMassHypertriton"), v0.mHypertriton());
            registry.fill(HIST("hMassHypertritonTotal"), v0.mHypertriton());
            registry.fill(HIST("h3dMassHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mHypertriton());
            registry.fill(HIST("h3dTotalHypertriton"), ct, v0.pt(), v0.mHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mHypertriton());
            }
          }
        }
      }

      // AntiHypertriton
      if (TMath::Abs( GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) ) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {

        registry.fill(HIST("hSelectedV0Counter"), 6.5);
        if(v0.posTrack_as<MyTracks>().pt() > 0.2 && v0.posTrack_as<MyTracks>().pt() < 1.2 && v0.negTrack_as<MyTracks>().pt() > 1.8 && v0.negTrack_as<MyTracks>().pt() < 10 && v0.pt() > 2 && v0.pt() < 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 7.5);
          if (TMath::Abs(v0.dcapostopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 8.5);

            registry.fill(HIST("hPtAntiHelium3"), v0.negativept());
            registry.fill(HIST("hPtPion"), v0.positivept());
            registry.fill(HIST("hMassAntiHypertriton"), v0.mAntiHypertriton());
            registry.fill(HIST("hMassHypertritonTotal"), v0.mAntiHypertriton());
            registry.fill(HIST("h3dMassAntiHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mAntiHypertriton());
            registry.fill(HIST("h3dTotalHypertriton"), ct, v0.pt(), v0.mAntiHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassAntiHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiHypertriton());
            }
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
      adaptAnalysisTask<hypertritonQa>(cfgc),
  };
}
