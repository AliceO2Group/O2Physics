// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
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
// Adaptation of V0 analysis task for run on MC data
// ========================
//
// This code loops over a V0Data table and produces some
// standard analysis output. It requires either
// the hypertritonfinder or the hypertritonbuilder tasks
// to have been executed in the workflow (before).
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
#include "Common/DataModel/PIDResponse.h"

#include "DetectorsVertexing/DCAFitterN.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>
#include "DetectorsBase/Propagator.h"
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

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,  aod::pidTPCFullPi, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::McTrackLabels>;

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
        {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
        {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
        {"hArmenterosPreAnalyserCuts", "hArmenterosPreAnalyserCuts", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}}},
        {"hMassHypertritonMCportion", "hMassHypertritonMCportion", {HistType::kTH1F, {{1000, 2.5f, 3.5f, "Inv. Mass (GeV/c^{2})"}}}},

      },
  };

  void init(InitContext const&)
  {
    AxisSpec massAxisHypertriton = {120, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {massAxisHypertriton}});
    registry.add("hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {massAxisHypertriton}});
  }

  void process(aod::Collision const& collision, aod::V0Datas const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
  {
    for (auto& v0 : fullV0s) {
      registry.fill(HIST("hMassHypertriton"), v0.mHypertriton());
      registry.fill(HIST("hMassAntiHypertriton"), v0.mAntiHypertriton());

      registry.fill(HIST("hV0Radius"), v0.v0radius());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
      registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());

      registry.fill(HIST("hArmenterosPreAnalyserCuts"), v0.alpha(), v0.qtarm());

      auto reconegtrack = v0.negTrack_as<MyTracks>();
      auto recopostrack = v0.posTrack_as<MyTracks>();
      if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle()) {
        continue;
      }
      auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
      auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();
      if (!mcnegtrack.has_mothers() || !mcpostrack.has_mothers()) {
        continue;
      }
      for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
        for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
          // LOGF(info, "MotherNeg = %d, MotherPos = %d", particleMotherOfNeg.pdgCode(), particleMotherOfPos.pdgCode());
          if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
            registry.fill(HIST("hMassHypertritonMCportion"), v0.mHypertriton());
          }
          if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
            registry.fill(HIST("hMassHypertritonMCportion"), v0.mAntiHypertriton());
          }
        }
      }

    }
  }
};

struct hypertritonAnalysisMc {

  HistogramRegistry registry{
    "registry",
      {
        {"h3dMassHypertriton", "h3dMassHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"h3dMassAntiHypertriton", "h3dMassAntiHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"h3dMassHypertriton_MC_truePt", "h3dMassHypertriton_MC_truePt", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"h3dMassAntiHypertriton_MC_truePt", "h3dMassAntiHypertriton_MC_truePt", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
        {"hMassHypertritonTotal", "hMassHypertritonTotal", {HistType::kTH1F, {{120, 2.9f, 3.2f}}}},

        {"hSelectedV0Counter", "hSelectedV0Counter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
        {"hTestCounter", "hTestCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},


        {"hSelectedEventCounter", "hSelectedEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},

        {"hArmenterosPostAnalyserCuts", "hArmenterosPostAnalyserCuts", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}}},
        {"hArmenterosPostAnalyserCuts_MC", "hArmenterosPostAnalyserCuts_MC", {HistType::kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}}}},
      },
  };
  //_MC_truePt histograms: additional .pdgcode cut, and fill true Pt instead of reconstructed Pt

  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};
  ConfigurableAxis massHypertritonbinning{"Hypertriton-mass-binning", {40, 2.95f, 3.05f}, ""};

  Configurable<int> saveDcaHist{"saveDcaHist", 1, "saveDcaHist"};
  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisHypertriton = {massHypertritonbinning, "Inv. Mass (GeV/c^{2})"};

    if (saveDcaHist == 1){
      registry.add("h3dMassHypertritonDca", "h3dMassHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
      registry.add("h3dMassAntiHypertritonDca", "h3dMassAntiHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
      registry.add("h3dMassHypertritonDca_MC_truePt", "h3dMassHypertritonDca_MC_truePt", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
      registry.add("h3dMassAntiHypertritonDca_MC_truePt", "h3dMassAntiHypertritonDca_MC_truePt", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
    }

    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(2, "V0CosPA");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(3, "hasMC");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(4, "TrackEta");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(5, "HasMC&Hypertriton Rapidity");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(6, "lifetime");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(7, "DcaV0Dau");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(8, "TPCPID");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(9, "PtCut");
    registry.get<TH1>(HIST("hSelectedV0Counter"))->GetXaxis()->SetBinLabel(10, "PionDcatoPV");
  }

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcapiontopv{"dcapiontopv", .1, "DCA Pion To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> etacut{"etacut", 0.9, "etacut"};
  Configurable<float> rapidity{"rapidity", 0.8, "rapidity"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  //Configurable<bool> boolArmenterosCut{"boolArmenterosCut", true, "cut on Armenteros-Podolanski graph"};
  //Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  //Configurable<bool> hasItsTest{"hasItsTest", false, "hasItsTest"};

  static constexpr float defaultLifetimeCuts[1][1] = {{40.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 1, {"lifetimecutHypertriton"}}, "lifetimecut"};

  //Filter preFilterV0 = aod::v0data::dcaV0daughters < dcav0dau;

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
    // void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
  {
    //Fix: the selection is not fine for MC data
    /*if (eventSelection && !collision.sel8()) {
      return;
      }*/
    registry.fill(HIST("hSelectedEventCounter"), 0.5);

    for (auto& v0 : fullV0s) {
      //   FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      registry.fill(HIST("hSelectedV0Counter"), 0.5);
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 1.5);
      if ( TMath::Abs(v0.posTrack_as<MyTracks>().eta()) > etacut || TMath::Abs(v0.negTrack_as<MyTracks>().eta()) > etacut){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 2.5);

      auto reconegtrack = v0.negTrack_as<MyTracks>();
      auto recopostrack = v0.posTrack_as<MyTracks>();
      if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle()) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 3.5);
      auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
      auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();

      if (TMath::Abs(v0.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 4.5);

      if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991 > lifetimecut->get("lifetimecutHypertriton")) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 5.5);
      if (v0.dcaV0daughters() > dcav0dau){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 6.5);

      //Hypertriton
      //if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
      if (TMath::Abs( GetTPCNSigmaHe3( 2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal()) ) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
        registry.fill(HIST("hSelectedV0Counter"), 7.5);

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
          registry.fill(HIST("hSelectedV0Counter"), 8.5);
          if (TMath::Abs(v0.dcanegtopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 9.5);

            registry.fill(HIST("hMassHypertritonTotal"), v0.mHypertriton());
            registry.fill(HIST("h3dMassHypertriton"), 0., v0.pt(), v0.mHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());

            for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
              for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
                if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
                  if (particleMotherOfNeg.isPhysicalPrimary()) {
                    registry.fill(HIST("h3dMassHypertriton_MC_truePt"), 0., particleMotherOfNeg.pt(), v0.mHypertriton());
                    registry.fill(HIST("hArmenterosPostAnalyserCuts_MC"), v0.alpha(), v0.qtarm());
                  }
                }
                if (saveDcaHist == 1) {
                  registry.fill(HIST("h3dMassHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mHypertriton());

                  if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
                    registry.fill(HIST("h3dMassHypertritonDca_MC_truePt"), v0.dcaV0daughters(), particleMotherOfNeg.pt(), v0.mHypertriton());
                  }
                }
              }
            }
          }
        }
      }

      // AntiHypertriton
      //if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
      if (TMath::Abs( GetTPCNSigmaHe3( 2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal()) )  < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {

        registry.fill(HIST("hSelectedV0Counter"), 7.5);

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
          registry.fill(HIST("hSelectedV0Counter"), 8.5);
          if (TMath::Abs(v0.dcapostopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 9.5);

            registry.fill(HIST("hMassHypertritonTotal"), v0.mAntiHypertriton());
            registry.fill(HIST("h3dMassAntiHypertriton"), 0., v0.pt(), v0.mAntiHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());

            for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
              for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
                if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
                  if (particleMotherOfNeg.isPhysicalPrimary()) {
                    registry.fill(HIST("h3dMassAntiHypertriton_MC_truePt"), 0., particleMotherOfNeg.pt(), v0.mAntiHypertriton());
                    registry.fill(HIST("hArmenterosPostAnalyserCuts_MC"), v0.alpha(), v0.qtarm());
                  }
                }
                if (saveDcaHist == 1) {
                  registry.fill(HIST("h3dMassAntiHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiHypertriton());

                  if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
                    registry.fill(HIST("h3dMassAntiHypertritonDca_MC_truePt"), v0.dcaV0daughters(), particleMotherOfNeg.pt(), v0.mAntiHypertriton());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(hypertritonAnalysisMc, processRun3, "Process Run 3 data", true);

  //Fix: Set as same as Run3!
  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
  {
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    registry.fill(HIST("hSelectedEventCounter"), 0.5);

    for (auto& v0 : fullV0s) {
      //   FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      registry.fill(HIST("hSelectedV0Counter"), 0.5);
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 1.5);
      if ( TMath::Abs(v0.posTrack_as<MyTracks>().eta()) > etacut || TMath::Abs(v0.negTrack_as<MyTracks>().eta()) > etacut){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 2.5);

      auto reconegtrack = v0.negTrack_as<MyTracks>();
      auto recopostrack = v0.posTrack_as<MyTracks>();
      if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle()) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 3.5);
      auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
      auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();

      if (TMath::Abs(v0.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 4.5);

      if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991 > lifetimecut->get("lifetimecutHypertriton")) {
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 5.5);
      if (v0.dcaV0daughters() > dcav0dau){
        continue;
      }
      registry.fill(HIST("hSelectedV0Counter"), 6.5);

      //Hypertriton
      //if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
      if (TMath::Abs( GetTPCNSigmaHe3( 2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal()) ) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
        registry.fill(HIST("hSelectedV0Counter"), 7.5);

        if(v0.negativept() > 0.2 && v0.negativept() < 1.2 && v0.positivept() > 1.8 && v0.positivept() < 10 && v0.pt() > 2 && v0.pt() < 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 8.5);
          if (TMath::Abs(v0.dcanegtopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 9.5);

            registry.fill(HIST("hMassHypertritonTotal"), v0.mHypertriton());
            registry.fill(HIST("h3dMassHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());

            for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
              for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
                if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
                  if (particleMotherOfNeg.isPhysicalPrimary()) {
                    registry.fill(HIST("h3dMassHypertriton_MC_truePt"), collision.centRun2V0M(), particleMotherOfNeg.pt(), v0.mHypertriton());
                    registry.fill(HIST("hArmenterosPostAnalyserCuts_MC"), v0.alpha(), v0.qtarm());
                  }
                }
                if (saveDcaHist == 1) {
                  registry.fill(HIST("h3dMassHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mHypertriton());

                  if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
                    registry.fill(HIST("h3dMassHypertritonDca_MC_truePt"), v0.dcaV0daughters(), particleMotherOfNeg.pt(), v0.mHypertriton());
                  }
                }
              }
            }
          }
        }
      }

      // AntiHypertriton
      //if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
      if (TMath::Abs( GetTPCNSigmaHe3( 2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal()) )  < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {

        registry.fill(HIST("hSelectedV0Counter"), 7.5);
        if(v0.positivept() > 0.2 && v0.positivept() < 1.2 && v0.negativept() > 1.8 && v0.negativept() < 10 && v0.pt() > 2 && v0.pt() < 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 8.5);
          if (TMath::Abs(v0.dcapostopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 9.5);

            registry.fill(HIST("hMassHypertritonTotal"), v0.mAntiHypertriton());
            registry.fill(HIST("h3dMassAntiHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mAntiHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());

            for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
              for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
                if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
                  if (particleMotherOfNeg.isPhysicalPrimary()) {
                    registry.fill(HIST("h3dMassAntiHypertriton_MC_truePt"), collision.centRun2V0M(), particleMotherOfNeg.pt(), v0.mAntiHypertriton());
                    registry.fill(HIST("hArmenterosPostAnalyserCuts_MC"), v0.alpha(), v0.qtarm());
                  }
                }
                if (saveDcaHist == 1) {
                  registry.fill(HIST("h3dMassAntiHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiHypertriton());

                  if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
                    registry.fill(HIST("h3dMassAntiHypertritonDca_MC_truePt"), v0.dcaV0daughters(), particleMotherOfNeg.pt(), v0.mAntiHypertriton());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(hypertritonAnalysisMc, processRun2, "Process Run 2 data", false);
};





namespace o2::aod
{
  namespace v0goodhe3
  {
    DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
    DECLARE_SOA_INDEX_COLUMN(Collision, collision);
  }
  DECLARE_SOA_TABLE(V0GoodHe3, "AOD", "V0GOODHE3", o2::soa::Index<>, v0goodhe3::GoodTrackId, v0goodhe3::CollisionId);
  namespace v0goodpion
  {
    DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
    DECLARE_SOA_INDEX_COLUMN(Collision, collision);
  } 
  DECLARE_SOA_TABLE(V0GoodPion, "AOD", "V0GOODPION", o2::soa::Index<>, v0goodpion::GoodTrackId, v0goodpion::CollisionId);
} 

struct hypertritonTrackMcinfo{
  //Basic checks
  HistogramRegistry registry{
    "registry",
      {

        {"hTotalCollCounter", "hTotalCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hParticleCount", "hParticleCount", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hParticleCount2", "hParticleCount2", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},//for tpcncls > 70
        {"hDauHelium3Count", "hDauHelium3Count", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hDauPionCount", "hDauPionCount", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}},
        {"hBackgroundPionCount", "hBackgroundPionCount", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
        {"hTPCNClsCrossedRows", "hTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
        {"hTrackEta", "hTrackEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hTrackITSNcls", "hTrackITSNcls", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
        {"hTrackMcRapidity", "hTrackMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hTrackNsigmaHelium3", "hTrackNsigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hTrackNsigmaPion", "hTrackNsigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHypertritonEta", "hHypertritomEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHypertritonMcRapidity", "hHypertritonMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHelium3Eta", "hHelium3Eta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHelium3McRapidity", "hHelium3McRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},


        {"hHypertritonMcPt", "hHypertritonMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},

        {"hHelium3Px", "hHelium3Px", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3Py", "hHelium3Py", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3Pz", "hHelium3Pz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3Pt", "hHelium3Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hHelium3P", "hHelium3P", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hHelium3McPx", "hHelium3McPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3McPy", "hHelium3McPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3McPz", "hHelium3McPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3McPt", "hHelium3McPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hHelium3McP", "hHelium3McP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hPionPx", "hPionPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionPy", "hPionPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionPz", "hPionPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionPt", "hPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hPionP", "hPionP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hPionMcPx", "hPionMcPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionMcPy", "hPionMcPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionMcPz", "hPionMcPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionMcPt", "hPionMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hPionMcP", "hPionMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},


        {"hHelium3ITSNcls", "hHelium3ITSNcls", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
        {"hHelium3NsigmaHelium3", "hHelium3NsigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3NsigmaPion", "hHelium3NsigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3NsigmaTriton", "hHelium3NsigmaTriton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3TPCBB", "hHelium3TPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal" }}}},
        {"hHelium3NsigmaHelium3Test", "hHelium3NsigmaHelium3Test", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3TPCBBTest", "hHelium3TPCBBTest", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal" }}}},
        {"hHelium3WrongNSigmaDCAToPV", "hHelium3WrongNSigmaDCAToPV", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHelium3WrongNSigmaPt", "hHelium3WrongNSigmaPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hHelium3WrongNSigmaP", "hHelium3WrongNSigmaP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},

        {"hHelium3TPCNClsCrossedRows", "hHelium3TPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
        {"hPionTPCNClsCrossedRows", "hPionTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},

        {"hTPCBB", "hTPCBB", {HistType::kTH2F, {{120, -8.0f, 8.0f, "p/z(GeV/c)"}, {100, 0.0f, 1000.0f, "TPCSignal" }}}},
      },
  };

  Produces<aod::V0GoodHe3> v0GoodHe3Tracks;
  Produces<aod::V0GoodPion> v0GoodPionTracks;

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(2, "Has_mcparticle");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(3, "Rapidity Cut(off)");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(4, "McisHypertriton");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(5, "McisHelium3");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(6, "McisPion");

    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(2, "Has_mcparticle");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(3, "Rapidity Cut(off)");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(4, "McisHypertriton");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(5, "McisHelium3");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(6, "McisPion");

    registry.get<TH1>(HIST("hDauHelium3Count"))->GetXaxis()->SetBinLabel(1, "hasMom");
    registry.get<TH1>(HIST("hDauHelium3Count"))->GetXaxis()->SetBinLabel(2, "FromHypertriton");
    registry.get<TH1>(HIST("hDauHelium3Count"))->GetXaxis()->SetBinLabel(3, "TPCNcls");
    registry.get<TH1>(HIST("hDauHelium3Count"))->GetXaxis()->SetBinLabel(4, "Eta");
    registry.get<TH1>(HIST("hDauHelium3Count"))->GetXaxis()->SetBinLabel(5, "Pt");
    registry.get<TH1>(HIST("hDauHelium3Count"))->GetXaxis()->SetBinLabel(6, "TPCPID");
    registry.get<TH1>(HIST("hDauPionCount"))->GetXaxis()->SetBinLabel(1, "hasMom");
    registry.get<TH1>(HIST("hDauPionCount"))->GetXaxis()->SetBinLabel(2, "FromHypertriton");
    registry.get<TH1>(HIST("hDauPionCount"))->GetXaxis()->SetBinLabel(3, "TPCNcls");
    registry.get<TH1>(HIST("hDauPionCount"))->GetXaxis()->SetBinLabel(4, "Eta");
    registry.get<TH1>(HIST("hDauPionCount"))->GetXaxis()->SetBinLabel(5, "Pt");
    registry.get<TH1>(HIST("hDauPionCount"))->GetXaxis()->SetBinLabel(6, "TPCPID");
    registry.get<TH1>(HIST("hDauPionCount"))->GetXaxis()->SetBinLabel(7, "DcatoPV");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(1, "NotFromHypertriton");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(2, "TPCNcls");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(3, "Eta");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(4, "Pt");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(5, "TPCPID");
  }


  void process(aod::Collision const& collision, aod::McParticles const& mcparticles, MyTracks const& tracks)
  {

    registry.fill(HIST("hTotalCollCounter"), 0.5);

    for (auto& track : tracks) {
      registry.fill(HIST("hParticleCount"), 0.5);
      registry.fill(HIST("hTrackITSNcls"), track.itsNCls());
      if (track.tpcNClsCrossedRows() > 70){
        registry.fill(HIST("hParticleCount2"), 0.5);
      }
      registry.fill(HIST("hTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("hTrackNsigmaHelium3"), GetTPCNSigmaHe3(2*track.p(), track.tpcSignal()) );
      registry.fill(HIST("hTrackNsigmaPion"), track.tpcNSigmaPi());

      if(!track.has_mcParticle()){
        continue;
      }
      registry.fill(HIST("hParticleCount"), 1.5);
      auto mcparticle = track.mcParticle_as<aod::McParticles>();
      if (track.tpcNClsCrossedRows() > 70){
        registry.fill(HIST("hParticleCount2"), 1.5);
      }
      registry.fill(HIST("hTPCBB"), track.p()*track.sign(), track.tpcSignal());

      //if (TMath::Abs(mcparticle.y()) > 0.9) {continue;}
      registry.fill(HIST("hParticleCount"), 2.5);
      if (track.tpcNClsCrossedRows() > 70){
        registry.fill(HIST("hParticleCount2"), 2.5);
      }
      registry.fill(HIST("hTrackEta"), track.eta());
      registry.fill(HIST("hTrackMcRapidity"), mcparticle.y());


      //Hypertriton
      if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hParticleCount"), 3.5);
        if (track.tpcNClsCrossedRows() > 70){
          registry.fill(HIST("hParticleCount2"), 3.5);
        }
        registry.fill(HIST("hHypertritonMcPt"), mcparticle.pt());
        registry.fill(HIST("hHypertritonEta"), track.eta());
        registry.fill(HIST("hHypertritonMcRapidity"), mcparticle.y());
      }

      //Helium3
      if (mcparticle.pdgCode() == 1000020030 || mcparticle.pdgCode() == -1000020030) {
        registry.fill(HIST("hParticleCount"), 4.5);
        if (track.tpcNClsCrossedRows() > 70){
          registry.fill(HIST("hParticleCount2"), 4.5);
          registry.fill(HIST("hHelium3NsigmaHelium3Test"), GetTPCNSigmaHe3(2*track.p(), track.tpcSignal()));
          registry.fill(HIST("hHelium3TPCBBTest"), track.p()*track.sign(), track.tpcSignal());
          if (TMath::Abs(GetTPCNSigmaHe3(2*track.p(), track.tpcSignal())) > 5){
            registry.fill(HIST("hHelium3WrongNSigmaDCAToPV"), track.dcaXY());
            registry.fill(HIST("hHelium3WrongNSigmaPt"), 2*track.pt());
            registry.fill(HIST("hHelium3WrongNSigmaP"), 2*track.p());
          }
        }

        if (mcparticle.has_mothers()){
          registry.fill(HIST("hDauHelium3Count"), 0.5);
          for (auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
            if ( particleMother.pdgCode() != 1010010030 && particleMother.pdgCode() != -1010010030){
              continue;
            }
            registry.fill(HIST("hDauHelium3Count"), 1.5);
            if (track.tpcNClsCrossedRows() < 70) {
              continue;
            }
            registry.fill(HIST("hDauHelium3Count"), 2.5);
            if (TMath::Abs(track.eta()) > 0.9) {
              continue;
            }
            registry.fill(HIST("hDauHelium3Count"), 3.5);
            if ( 2*track.pt() < 1.8 || 2*track.pt() > 10) {
              continue;
            }
            registry.fill(HIST("hDauHelium3Count"), 4.5);
            if (TMath::Abs(GetTPCNSigmaHe3(2*track.p(), track.tpcSignal()) ) > 5) {
              continue;
            }
            registry.fill(HIST("hDauHelium3Count"), 5.5);
            registry.fill(HIST("hHelium3ITSNcls"), track.itsNCls());
            v0GoodHe3Tracks(track.globalIndex(), track.collisionId());
          }
        }

        registry.fill(HIST("hHelium3McPx"), mcparticle.px());
        registry.fill(HIST("hHelium3McPy"), mcparticle.py());
        registry.fill(HIST("hHelium3McPz"), mcparticle.pz());
        registry.fill(HIST("hHelium3McPt"), mcparticle.pt());
        registry.fill(HIST("hHelium3McP"), mcparticle.p());

        registry.fill(HIST("hHelium3Px"), 2*track.px());
        registry.fill(HIST("hHelium3Py"), 2*track.py());
        registry.fill(HIST("hHelium3Pz"), 2*track.pz());
        registry.fill(HIST("hHelium3Pt"), 2*track.pt());
        registry.fill(HIST("hHelium3P"), 2*track.p());

        registry.fill(HIST("hHelium3NsigmaHelium3"), GetTPCNSigmaHe3(2*track.p(), track.tpcSignal()));
        registry.fill(HIST("hHelium3NsigmaPion"), track.tpcNSigmaPi());
        registry.fill(HIST("hHelium3NsigmaTriton"), track.tpcNSigmaTr());
        registry.fill(HIST("hHelium3TPCNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hHelium3Eta"), track.eta());
        registry.fill(HIST("hHelium3McRapidity"), mcparticle.y());
        registry.fill(HIST("hHelium3TPCBB"), track.p()*track.sign(), track.tpcSignal());

      }

      //Pion
      if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
        registry.fill(HIST("hParticleCount"), 5.5);
        if (track.tpcNClsCrossedRows() > 70){
          registry.fill(HIST("hParticleCount2"), 5.5);
        }

        if (mcparticle.has_mothers()){
          registry.fill(HIST("hDauPionCount"), 0.5);
          for (auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
            if ( particleMother.pdgCode() != 1010010030 && particleMother.pdgCode() != -1010010030){
              continue;
            }
            registry.fill(HIST("hDauPionCount"), 1.5);
            if (track.tpcNClsCrossedRows() < 70) {
              continue;
            }
            registry.fill(HIST("hDauPionCount"), 2.5);
            if (TMath::Abs(track.eta()) > 0.9) {
              continue;
            }
            registry.fill(HIST("hDauPionCount"), 3.5);
            if ( track.pt() < 0.2 || track.pt() > 1.2) {
              continue;
            }
            registry.fill(HIST("hDauPionCount"), 4.5);
            if (TMath::Abs(track.tpcNSigmaPi()) > 5) {
              continue;
            }
            registry.fill(HIST("hDauPionCount"), 5.5);
            if ( TMath::Abs(track.dcaXY()) < 0.1) {
              continue;
            }
            registry.fill(HIST("hDauPionCount"), 6.5);
            v0GoodPionTracks(track.globalIndex(), track.collisionId());
          }
        }

        bool isFromHypertriton = false; 
        if (mcparticle.has_mothers()){
          for (auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
            if ( particleMother.pdgCode() == 1010010030 || particleMother.pdgCode() == -1010010030){
              isFromHypertriton = true;
            }
          }
        }

        while (!isFromHypertriton){
            registry.fill(HIST("hBackgroundPionCount"), 0.5);
            if (track.tpcNClsCrossedRows() < 70) {
              break;
            }
            registry.fill(HIST("hBackgroundPionCount"), 1.5);
            if (TMath::Abs(track.eta()) > 0.9) {
              break;
            }
            registry.fill(HIST("hBackgroundPionCount"), 2.5);
            if ( track.pt() < 0.2 || track.pt() > 1.2) {
              break;
            }
            registry.fill(HIST("hBackgroundPionCount"), 3.5);
            if (TMath::Abs(track.tpcNSigmaPi()) > 5) {
              break;
            }
            registry.fill(HIST("hBackgroundPionCount"), 4.5);
              break;
        }

        registry.fill(HIST("hPionMcPx"), mcparticle.px());
        registry.fill(HIST("hPionMcPy"), mcparticle.py());
        registry.fill(HIST("hPionMcPz"), mcparticle.pz());
        registry.fill(HIST("hPionMcPt"), mcparticle.pt());
        registry.fill(HIST("hPionMcP"), mcparticle.p());

        registry.fill(HIST("hPionPx"), track.px());
        registry.fill(HIST("hPionPy"), track.py());
        registry.fill(HIST("hPionPz"), track.pz());
        registry.fill(HIST("hPionPt"), track.pt());
        registry.fill(HIST("hPionP"), track.p());
        registry.fill(HIST("hPionTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
      }
    }
  }
};

struct goodtrackcheck{
  // Configurables
  Configurable<double> d_UseAbsDCA{"d_UseAbsDCA", kTRUE, "Use Abs DCAs"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};

  HistogramRegistry registry{
    "registry",
      {
        {"hV0CutCounter", "hV0CutCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
        {"hTrueV0Counter", "hTrueV0Counter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
        {"hTrueV0PhysicalCounter", "hTrueV0PhysicalCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
        {"hHe3Counter", "hHe3Counter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hPionCounter", "hPionCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hMcMassHypertriton", "hMcMassHypertriton", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
      },
  };

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  int mRunNumber;
  float d_bz;
  float maxSnp;  //max sine phi for propagation
  float maxStep; //max step size (cm) for propagation
  void init(InitContext& context)
  {
    // using namespace analysis::lambdakzerobuilder;
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  //could be changed later
    maxStep = 2.00f; //could be changed later

    ccdb->setURL("https://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/Geometry");
      /* it seems this is needed at this level for the material LUT to work properly */
      /* but what happens if the run changes while doing the processing?             */
      constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(1, "Sign");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(2, "DiffCol");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(3, "hasSV");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(4, "hasSV2");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(5, "Dcav0Dau");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(6, "CosPA");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(7, "Rapidity");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(8, "Lifetime");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(9, "pT");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(1, "Sign");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(2, "DiffCol");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(3, "hasSV");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(4, "hasSV2");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(5, "Dcav0Dau");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(6, "CosPA");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(7, "Rapidity");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(8, "Lifetime");
    registry.get<TH1>(HIST("hTrueV0Counter"))->GetXaxis()->SetBinLabel(9, "pT");
  }

  float getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    float output = grpo->getNominalL3Field();
    return output;
  }

  void CheckAndUpdate(Int_t lRunNumber, uint64_t lTimeStamp)
  {
    if (lRunNumber != mRunNumber) {
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(lTimeStamp);
      } else {
        d_bz = d_bz_input;
      }
      mRunNumber = lRunNumber;
    }
  }

  void process(aod::Collision const& collision, MyTracks const& tracks, aod::McParticles const& mcparticles, 
      aod::V0GoodHe3 const& he3tracks, aod::V0GoodPion const& piontracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    CheckAndUpdate(bc.runNumber(), bc.timestamp());

    // Define o2 fitter, 2-prong
    o2::vertexing::DCAFitterN<2> fitter;
    fitter.setBz(d_bz);
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);


    for (auto& t0id : he3tracks) {
      registry.fill(HIST("hHe3Counter"), 0.5);
    }
    for (auto& t1id : piontracks) {
      registry.fill(HIST("hPionCounter"), 0.5);
    }

    for (auto& t0id : he3tracks) { // FIXME: turn into combination(...)
      for (auto& t1id : piontracks) {


        auto t0 = t0id.goodTrack_as<MyTracks>();
        auto t1 = t1id.goodTrack_as<MyTracks>();
        auto t0mc = t0.mcParticle_as<aod::McParticles>();
        auto t1mc = t1.mcParticle_as<aod::McParticles>();
        auto Track1 = getTrackParCov(t0);
        auto Track2 = getTrackParCov(t1);
        auto he3Track = getTrackParCov(t0);
        auto pionTrack = getTrackParCov(t1);

        if (t0.sign() + t1.sign() != 0.0f){
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 0.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 0.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 0.5);
              }
            }
          }
        }
        if (t0.collisionId() != t1.collisionId()) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 1.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 1.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 1.5);
              }
            }
          }
        }

        // Try to progate to dca
        int nCand = fitter.process(Track1, Track2);
        if (nCand == 0) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 2.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 2.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 2.5);
              }
            }
          }
        }

        //------------------copy from lamdakzerobuilder---------------------
        double finalXpos = fitter.getTrack(0).getX();
        double finalXneg = fitter.getTrack(1).getX();

        // Rotate to desired alpha
        he3Track.rotateParam(fitter.getTrack(0).getAlpha());
        pionTrack.rotateParam(fitter.getTrack(1).getAlpha());

        // Retry closer to minimum with material corrections
        o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
        if (useMatCorrType == 1)
          matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
        if (useMatCorrType == 2)
          matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

        o2::base::Propagator::Instance()->propagateToX(he3Track, finalXpos, d_bz, maxSnp, maxStep, matCorr);
        o2::base::Propagator::Instance()->propagateToX(pionTrack, finalXneg, d_bz, maxSnp, maxStep, matCorr);

        nCand = fitter.process(he3Track, pionTrack);
        if (nCand == 0) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 3.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 3.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 3.5);
              }
            }
          }
        }

        //------------------------------------------------------------------

        const auto& vtx = fitter.getPCACandidate();

        // DCA V0 daughters
        auto thisdcav0dau = fitter.getChi2AtPCACandidate();
        if (thisdcav0dau > dcav0dau) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 4.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 4.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 4.5);
              }
            }
          }
        }

        std::array<float, 3> pos = {0.};
        std::array<float, 3> pvec0;
        std::array<float, 3> pvec1;
        for (int i = 0; i < 3; i++) {
          pos[i] = vtx[i];
        }
        //fitter.getTrack(0).getPxPyPzGlo(pvec0);
        //fitter.getTrack(1).getPxPyPzGlo(pvec1);

        //------------------copy from lamdakzerobuilder---------------------
        he3Track.getPxPyPzGlo(pvec0);
        pionTrack.getPxPyPzGlo(pvec1);
        //------------------------------------------------------------------
        int pTrackCharge = 1, nTrackCharge = 1;
      if (TMath::Abs( GetTPCNSigmaHe3( 2*t0.p(), t0.tpcSignal()) ) < 5){
          pTrackCharge = 2;
        } 
      if (TMath::Abs( GetTPCNSigmaHe3( 2*t1.p(), t1.tpcSignal()) ) < 5){
          nTrackCharge = 2;
        } 
        for (int i=0; i<3; i++){
          pvec0[i] = pvec0[i] * pTrackCharge;
          pvec1[i] = pvec1[i] * nTrackCharge;
        }

        auto thisv0cospa = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()},
            array{vtx[0], vtx[1], vtx[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
        if (thisv0cospa < v0cospa) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 5.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 5.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 5.5);
              }
            }
          }
        }

        double hypertritonMcMass = RecoDecay::m(array{array{t0mc.px(), t0mc.py(), t0mc.pz()}, array{t1mc.px(), t1mc.py(), t1mc.pz()}}, array{2.80839, RecoDecay::getMassPDG(kPiPlus)}); 
        double hypertritonMass = RecoDecay::m(array{array{pvec0[0], pvec0[1], pvec0[2]}, array{pvec1[0], pvec1[1], pvec1[2]}}, array{2.80839, RecoDecay::getMassPDG(kPiPlus)}); 
        double hypertritonRapidity = RecoDecay::y(array{pvec0[0]+pvec1[0], pvec0[1]+pvec1[1], pvec0[2]+pvec1[2]}, 2.991); 
        double hypertritonPt = TMath::Sqrt(TMath::Power(pvec0[0]+pvec1[0], 2) + TMath::Power(pvec0[1]+pvec1[1], 2));
        double ct = std::sqrt(std::pow(collision.posX() - pos[0], 2) + std::pow(collision.posY() - pos[1], 2) + std::pow(collision.posZ() - pos[2], 2)) / (TMath::Sqrt( TMath::Power(pvec0[0]+pvec1[0], 2) + TMath::Power( pvec0[1]+pvec1[1], 2) + TMath::Power( pvec0[2]+pvec1[2], 2)) + 1E-10) * 2.991;

        if (TMath::Abs(hypertritonRapidity) > 0.8){
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 6.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 6.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 6.5);
              }
            }
          }
        }

        if ( ct > 40){
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 7.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 7.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 7.5);
              }
            }
          }
        }

        if (hypertritonPt < 2 || hypertritonPt > 9){
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 8.5);
        for (auto& particleMother1 : t0mc.mothers_as<aod::McParticles>()) {
          for (auto& particleMother2 : t1mc.mothers_as<aod::McParticles>()) {
            if (particleMother1 == particleMother2 && TMath::Abs(particleMother1.pdgCode()) == 1010010030) {
              registry.fill(HIST("hTrueV0Counter"), 8.5);
              if (particleMother1.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueV0PhysicalCounter"), 8.5);
              }
            }
          }
        }

        registry.fill(HIST("hMassHypertriton"), hypertritonMass);
        registry.fill(HIST("hMcMassHypertriton"), hypertritonMcMass);
      }
    }
  }

};


namespace o2::aod
{
  namespace mcparticleidcheck{
    DECLARE_SOA_COLUMN(HypertritonDaughterPdgCode, hypertritonDaughterPdgCode, int);
    DECLARE_SOA_COLUMN(Helium3DaughterPdgCode, helium3DaughterPdgCode, int);
    DECLARE_SOA_COLUMN(Helium3MontherPdgCode, helium3MontherPdgCode, int);
  }
  DECLARE_SOA_TABLE(McHypertritonCheck, "AOD", "HypCheck", 
      mcparticleidcheck::HypertritonDaughterPdgCode 
      )
    DECLARE_SOA_TABLE(McHelium3Check, "AOD", "Helium3Check", 
        mcparticleidcheck::Helium3MontherPdgCode,
        mcparticleidcheck::Helium3DaughterPdgCode
        )
}


struct hypertritonMcParticleCount {
  //Basic checks
  Produces<aod::McHypertritonCheck> hypertrtitonCheckTable;
  Produces<aod::McHelium3Check> helium3CheckTable;
  HistogramRegistry registry{
    "registry",
      {
        {"hTotalMcCollCounter", "hTotalMcCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hSelAndRecoMcCollCounter", "hSelAndRecoMcCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},

        {"h3dMCHypertriton", "h3dMCHypertriton", {HistType::kTH3F, {{20, -1.0f, 1.0f, "Rapidity"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {50, 0.0f, 50.0f, "ct(cm)"}}}},
        {"hMcHypertritonCheck", "hMcHypertritonCheck", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
        {"hMcParticleCount", "hMcParticleCount", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hHypertritonCount_PtDiff", "hHypertritonCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},


        {"hMcHypertritonCount", "hMcHypertritonCount", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
        {"hMcHypertritonPt", "hMcHypertritonPt", {HistType::kTH1F, {{300, 0.0f, 15.0f}}}},
        {"hMcHypertritonLifeTime", "hMcHypertritonLifeTime", {HistType::kTH1F, {{500, 0.0f, 50.0f}}}},

        {"hMcHelium3Check", "hMcHelium3Check", {HistType::kTH1F, {{1, 0.0f, 4.0f}}}},
        {"hMcHelium3Pt", "hMcHelium3Pt", {HistType::kTH1F, {{300, 0.0f, 15.0f}}}},
        {"hMcPionPt", "hMcPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},

        {"hMcDauHelium3Vx", "hMcDauHelium3Vx", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hMcDauHelium3Vy", "hMcDauHelium3Vy", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hMcDauHelium3Vz", "hMcDauHelium3Vz", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(1, "PdgCode");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(2, "Rapidity");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(3, "Lifetime");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(4, "PtCut");

    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(1, "ReadinAfterEvtRecSelCut(off)");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(2, "y<0.9(off)");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(3, "IsPhysicalPrimary");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(4, "(Anti)Helium3");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(5, "(Anti)Hypertriton");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(6, "HasDaughter");
    registry.get<TH1>(HIST("hMcHypertritonCount"))->GetXaxis()->SetBinLabel(1, "Generated number");
    registry.get<TH1>(HIST("hMcHypertritonCount"))->GetXaxis()->SetBinLabel(2, "primary Hypertriton mothers");
    registry.get<TH1>(HIST("hMcHypertritonCount"))->GetXaxis()->SetBinLabel(3, "decaying into V0");
    registry.get<TH1>(HIST("hMcHelium3Check"))->GetXaxis()->SetBinLabel(1, "PdgCode");
  }

  Configurable<float> rapidityMCcut{"rapidityMCcut", 0.8, "rapidity cut MC count"};
  Configurable<bool> eventSelectionMC{"eventSelectionMC", true, "event selection MC count"};

  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>& collisions)
  {
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      /*if (eventSelectionMC && !collision.sel8()) {
        continue;
        }*/
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);

    registry.fill(HIST("hTotalMcCollCounter"), 0.5);

    /*const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();
      if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
      }*/
    registry.fill(HIST("hSelAndRecoMcCollCounter"), 0.5);

    for (auto& mcparticle : mcParticles) {

      if (mcparticle.pdgCode() == 1000020030 || mcparticle.pdgCode() == -1000020030) {
        registry.fill(HIST("hMcHelium3Check"), 0.5);
        registry.fill(HIST("hMcHelium3Pt"), mcparticle.pt());
        for (auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          for (auto& mcparticleMonther : mcparticle.mothers_as<aod::McParticles>()) {
            helium3CheckTable(mcparticleMonther.pdgCode(),mcparticleDaughter.pdgCode());
          }
        }
      }

      if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
        registry.fill(HIST("hMcPionPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hMcHypertritonCount"), 0.5);
        registry.fill(HIST("hMcHypertritonCheck"), 0.5);
        registry.fill(HIST("hMcHypertritonPt"), mcparticle.pt());

        double dauHe3Vx = -999; 
        double dauHe3Vy = -999;
        double dauHe3Vz = -999;
        for (auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          hypertrtitonCheckTable(mcparticleDaughter.pdgCode());
          if (mcparticleDaughter.pdgCode() == 1000020030 || mcparticleDaughter.pdgCode() == -1000020030){
            dauHe3Vx = mcparticleDaughter.vx();
            dauHe3Vy = mcparticleDaughter.vy();
            dauHe3Vz = mcparticleDaughter.vz();
          }
        }
        registry.fill(HIST("hMcDauHelium3Vx"), dauHe3Vx);
        registry.fill(HIST("hMcDauHelium3Vy"), dauHe3Vy);
        registry.fill(HIST("hMcDauHelium3Vz"), dauHe3Vz);

        double MClifetime = RecoDecay::sqrtSumOfSquares(dauHe3Vx - mcparticle.vx(), dauHe3Vy - mcparticle.vy(), dauHe3Vz - mcparticle.vz())*2.991/mcparticle.p();  
        registry.fill(HIST("hMcHypertritonLifeTime"), MClifetime);
        registry.fill(HIST("h3dMCHypertriton"), mcparticle.y(), mcparticle.pt(), MClifetime);

        //Count for hypertriton N_gen
        if (TMath::Abs(mcparticle.y()) < 0.8){
          registry.fill(HIST("hMcHypertritonCheck"), 1.5);
          if (MClifetime < 40){
            registry.fill(HIST("hMcHypertritonCheck"), 2.5);
            if (mcparticle.pt() > 2 && mcparticle.pt() < 9){
              registry.fill(HIST("hMcHypertritonCheck"), 3.5);
            }
          }
        }
      }

      registry.fill(HIST("hMcParticleCount"), 0.5);
      //if (TMath::Abs(mcparticle.y()) > rapidityMCcut) {continue;}
      registry.fill(HIST("hMcParticleCount"), 1.5);
      if (!mcparticle.isPhysicalPrimary()) {
        continue;
      }
      registry.fill(HIST("hMcParticleCount"), 2.5);
      if (mcparticle.pdgCode() == 1000020030 || mcparticle.pdgCode() == -1000020030) {
        registry.fill(HIST("hMcParticleCount"), 3.5);
      }
      if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hMcParticleCount"), 4.5);
      }
      if (!mcparticle.has_daughters()) {
        continue;
      }
      registry.fill(HIST("hMcParticleCount"), 5.5);

      if (std::abs(mcparticle.pdgCode()) == 1010010030) {
        registry.fill(HIST("hMcHypertritonCount"), 1.5);
        for (auto& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
          for (auto& mcparticleDaughter1 : mcparticle.daughters_as<aod::McParticles>()) {
            if ( (mcparticleDaughter0.pdgCode() == -211 && mcparticleDaughter1.pdgCode() == 1000020030) || (mcparticleDaughter0.pdgCode() == 211 && mcparticleDaughter1.pdgCode() == -1000020030) ) {
              registry.fill(HIST("hMcHypertritonCount"), 2.5);
              registry.fill(HIST("hHypertritonCount_PtDiff"), mcparticle.pt());
            }
          }
        }
      }
    }
  }
};

namespace o2::aod
{
  namespace v0mc{
    DECLARE_SOA_COLUMN(PosTrackMcPdgCode, posTrackMcPdgCode, int);
    DECLARE_SOA_COLUMN(PosTrackPdgPID, posTrackPdgPID, int);
    DECLARE_SOA_COLUMN(PosTrackPx, posTrackPx, float);
    DECLARE_SOA_COLUMN(PosTrackPy, posTrackPy, float);
    DECLARE_SOA_COLUMN(PosTrackPz, posTrackPz, float);
    DECLARE_SOA_COLUMN(PosTrackMcPx, posTrackMcPx, float);
    DECLARE_SOA_COLUMN(PosTrackMcPy, posTrackMcPy, float);
    DECLARE_SOA_COLUMN(PosTrackMcPz, posTrackMcPz, float);
    DECLARE_SOA_COLUMN(PosTrackTPCNCls, posTrackTPCNCls, int);
    DECLARE_SOA_COLUMN(PosTrackTPCSignals, posTrackTPCSignals, float);
    DECLARE_SOA_COLUMN(PosTrackTPCNSigmaHelium3, posTrackTPCNSigmaHe, float);
    DECLARE_SOA_COLUMN(PosTrackTPCNSigmaPion, posTrackTPCNSigmaPi, float);
    DECLARE_SOA_COLUMN(PosTrackTPCNSigmaTriton, posTrackTPCNSigmaTr, float);
    DECLARE_SOA_COLUMN(PosTrackDcaXY, posTrackDcaXY, float);

    DECLARE_SOA_COLUMN(NegTrackMcPdgCode, negTrackMcPdgCode, int);
    DECLARE_SOA_COLUMN(NegTrackPdgPID, negTrackPdgPID, int);
    DECLARE_SOA_COLUMN(NegTrackPx, negTrackPx, float);
    DECLARE_SOA_COLUMN(NegTrackPy, negTrackPy, float);
    DECLARE_SOA_COLUMN(NegTrackPz, negTrackPz, float);
    DECLARE_SOA_COLUMN(NegTrackMcPx, negTrackMcPx, float);
    DECLARE_SOA_COLUMN(NegTrackMcPy, negTrackMcPy, float);
    DECLARE_SOA_COLUMN(NegTrackMcPz, negTrackMcPz, float);
    DECLARE_SOA_COLUMN(NegTrackTPCNCls, negTrackTPCNCls, int);
    DECLARE_SOA_COLUMN(NegTrackTPCSignals, negTrackTPCSignals, int);
    DECLARE_SOA_COLUMN(NegTrackTPCNSigmaHelium3, negTrackTPCNSigmaHe, float);
    DECLARE_SOA_COLUMN(NegTrackTPCNSigmaPion, negTrackTPCNSigmaPion, float);
    DECLARE_SOA_COLUMN(NegTrackTPCNSigmaTriton, negTrackTPCNSigmaTriton, float);
    DECLARE_SOA_COLUMN(NegTrackDcaXY, negTrackDcaXY, float);
  }
  DECLARE_SOA_TABLE(V0Mc, "AOD", "v0mc", 
      v0mc::PosTrackMcPdgCode, 
      v0mc::PosTrackPdgPID, 
      v0mc::PosTrackPx, 
      v0mc::PosTrackPy, 
      v0mc::PosTrackPz, 
      v0mc::PosTrackMcPx,
      v0mc::PosTrackMcPy,
      v0mc::PosTrackMcPz,
      v0mc::PosTrackTPCNCls,
      v0mc::PosTrackTPCSignals,
      v0mc::PosTrackTPCNSigmaHelium3, 
      v0mc::PosTrackTPCNSigmaPion, 
      v0mc::PosTrackTPCNSigmaTriton, 
      v0mc::PosTrackDcaXY, 
      v0mc::NegTrackMcPdgCode, 
      v0mc::NegTrackPdgPID, 
      v0mc::NegTrackPx, 
      v0mc::NegTrackPy, 
      v0mc::NegTrackPz, 
      v0mc::NegTrackMcPx, 
      v0mc::NegTrackMcPy, 
      v0mc::NegTrackMcPz, 
      v0mc::NegTrackTPCNCls,
      v0mc::NegTrackTPCSignals,
      v0mc::NegTrackTPCNSigmaHelium3, 
      v0mc::NegTrackTPCNSigmaPion, 
      v0mc::NegTrackTPCNSigmaTriton, 
      v0mc::NegTrackDcaXY
        );

}

struct V0McCheck {

  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};//loose cut
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<float> lifetimecut{"lifetimecut", 40., "lifetimecut"}; //ct*/
  Configurable<float> etacut{"etacut", 0.9, "etacut"};
  Configurable<float> rapidity{"rapidity", 0.8, "rapidity"};
  Configurable<float> dcapiontopv{"dcapiontopv", .1, "DCA Pion To PV"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  //max sine phi for propagation
  float maxStep; //max step size (cm) for propagation

  Produces<aod::V0Mc> v0mcinfo;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{
    "registry",
      {
        {"hV0McCheckCounter", "hV0McCheckCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
        {"hTrueHypertritonCounter", "hTrueHypertritonCounter", {HistType::kTH1F, {{12, 0.0f, 12.0f}}}},
        {"hV0Helium3PerEvent", "hV0Helium3PerEvent", {HistType::kTH1F, {{7, -0.5f, 6.5f}}}},
        {"hV0AntiHelium3PerEvent", "hV0AntiHelium3PerEvent", {HistType::kTH1F, {{7, -0.5f, 6.5f}}}},
        {"hV0CollisionID", "hV0CollisionID", {HistType::kTH1F, {{100, 0.0f, 100.0f}}}},
        {"hV0PostrackPt", "hV0PostrackPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0NegtrackPt", "hV0PostrackPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0PostrackTPCNClus", "hV0PostrackTPCNClus", {HistType::kTH1F, {{24, 0.0f, 240.0f}}}},
        {"hV0NegtrackTPCNClus", "hV0NegtrackTPCNClus", {HistType::kTH1F, {{24, 0.0f, 240.0f}}}},
        {"hV0PostrackDcatoPv", "hV0PostrackDcatoPv", {HistType::kTH1F, {{10, 0.0f, 1.0f}}}},
        {"hV0NegtrackDcatoPv", "hV0NegtrackDcatoPv", {HistType::kTH1F, {{10, 0.0f, 1.0f}}}},
        {"hV0McPostrackPt", "hV0McPostrackPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0McNegtrackPt", "hV0McNegtrackPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0PostrackCharge", "hV0PostrackCharge", {HistType::kTH1F, {{10, -5.0f, 5.0f}}}},
        {"hV0NegtrackCharge", "hV0NegtrackCharge", {HistType::kTH1F, {{10, -5.0f, 5.0f}}}},
        {"hV0PostrackID", "hV0PostrackID", {HistType::kTH1F, {{20, 0.0f, 20.0f}}}},
        {"hV0NegtrackID", "hV0NegtrackID", {HistType::kTH1F, {{20, 0.0f, 20.0f}}}},

        {"hNSigmaHelium3Pos", "hNSigmaHelium3Pos", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaPionPos", "hNSigmaPionPos", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaTritonPos", "hNSigmaTritonPos", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaHelium3Neg", "hNSigmaHelium3Neg", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaPionNeg", "hNSigmaPionNeg", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hNSigmaTritonNeg", "hNSigmaTritonNeg", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hV0DausTPCBB", "hV0DausTPCBB", {HistType::kTH2F, {{120, -8.0f, 8.0f, "p/z(GeV/c)"}, {100, 0.0f, 1000.0f, "TPCSignal" }}}},

        {"hV0DauHelium3Px", "hV0DauHelium3Px", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3Py", "hV0DauHelium3Py", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3Pz", "hV0DauHelium3Pz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3Pt", "hV0DauHelium3Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0DauHelium3P", "hV0DauHelium3P", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0DauHelium3McPx", "hV0DauHelium3McPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3McPy", "hV0DauHelium3McPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3McPz", "hV0DauHelium3McPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3McPt", "hV0DauHelium3McPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0DauHelium3McP", "hV0DauHelium3McP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPx", "hV0DauHelium3DiffPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPy", "hV0DauHelium3DiffPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPz", "hV0DauHelium3DiffPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPt", "hV0DauHelium3DiffPt", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffP", "hV0DauHelium3DiffP", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3TPCBB", "hV0DauHelium3TPCBB", {HistType::kTH2F, {{120, -8.0f, 8.0f, "p/z(GeV/c)"}, {100, 0.0f, 1000.0f, "TPCSignal" }}}},
        /*{"hV0PionPx", "hV0PionPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionPy", "hV0PionPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionPz", "hV0PionPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionPt", "hV0PionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
          {"hV0PionP", "hV0PionP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
          {"hV0PionMcPx", "hV0PionMcPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionMcPy", "hV0PionMcPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionMcPz", "hV0PionMcPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionMcPt", "hV0PionMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
          {"hV0PionMcP", "hV0PionMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},*/
        {"hV0DauHelium3NSigmaHelium3", "hV0DauHelium3NSigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hV0DauHelium3NSigmaPion", "hV0DauHelium3NSigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hV0DauHelium3NSigmaTriton", "hV0DauHelium3NSigmaTriton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hV0DauPionNSigmaHelium3", "hV0DauPionNSigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hV0DauPionNSigmaPion", "hV0DauPionNSigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hV0McHypertritonMass", "hV0McHypertritonMass", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass", "hV0HypertritonMass", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass1", "hV0HypertritonMass1", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass2", "hV0HypertritonMass2", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass3", "hV0HypertritonMass3", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0McHypertritonMassAfterV0Cut", "hV0McHypertritonMassAfterV0Cut", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass4", "hV0HypertritonMass4", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass5", "hV0HypertritonMass5", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass6", "hV0HypertritonMass6", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass7", "hV0HypertritonMass7", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass8", "hV0HypertritonMass8", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass9", "hV0HypertritonMass9", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass10", "hV0HypertritonMass10", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass11", "hV0HypertritonMass11", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMass12", "hV0HypertritonMass12", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},
        {"hV0HypertritonMassTest", "hV0HypertritonMassTest", {HistType::kTH1F, {{80, 2.90f, 3.1f}}}},//for perfomance of v0radius
      },
  };

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  //could be changed later
    maxStep = 2.00f; //could be changed later

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
      /* it seems this is needed at this level for the material LUT to work properly */
      /* but what happens if the run changes while doing the processing?             */
      constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }

    registry.get<TH1>(HIST("hV0McCheckCounter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hV0McCheckCounter"))->GetXaxis()->SetBinLabel(2, "DauhasMc");

    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(1, "HypertritonPreCut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(2, "PhysicalHypertritonPreCut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(3, "HypertritonV0Cut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(4, "PhysicalHypertritonV0Cut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(5, "HypertritonCandidateCut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(6, "PhysicalHypertritonCandidateCut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(7, "AntiHypertritonPreCut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(8, "PhysicalAntiHypertritonPreCut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(9, "AntiHypertritonV0Cut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(10, "PhysicalAntiHypertritonV0Cut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(11, "AntiHypertritonCandidateCut");
    registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(12, "PhysicalAntiHypertritonCandidateCut");
  }

  float getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    float output = grpo->getNominalL3Field();
    return output;
  }

  void CheckAndUpdate(Int_t lRunNumber, uint64_t lTimeStamp)
  {
    if (lRunNumber != mRunNumber) {
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(lTimeStamp);
      } else {
        d_bz = d_bz_input;
      }
      mRunNumber = lRunNumber;
    }
  }
  void process(aod::Collision const& collision, aod::V0s const& V0s, MyTracks const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& mcparticles){

    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    CheckAndUpdate(bc.runNumber(), bc.timestamp());

    o2::vertexing::DCAFitterN<2> fitter;
    fitter.setBz(d_bz);
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true); // use d_UseAbsDCA once we want to use the weighted DCA

    int countHe3=0, countAntiHe3=0;
    std::vector<int> He3TrackID;
    std::vector<int> AntiHe3TrackID;
    for ( auto& v0 : V0s ){
      auto recopostrack = v0.posTrack_as<MyTracks>();
      auto reconegtrack = v0.negTrack_as<MyTracks>();
      if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle()) {
        continue;
      }
      auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();
      auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
      if (mcpostrack.pdgCode() == 1000020030 ){
        bool isNewHe3=true;;
        for( int id : He3TrackID){
          if( v0.posTrackId() == id){
            isNewHe3 = false;
            break;
          }
        }
        if (isNewHe3){
          countHe3++;
          He3TrackID.push_back(v0.posTrackId());
        }
      }

      if (mcnegtrack.pdgCode() == -1000020030 ){
        bool isNewAntiHe3=true;;
        for( int id : AntiHe3TrackID){
          if( v0.negTrackId() == id){
            isNewAntiHe3 = false;
            break;
          }
        }
        if (isNewAntiHe3){
          countAntiHe3++;
          AntiHe3TrackID.push_back(v0.negTrackId());
        }
      }

    }
    registry.fill(HIST("hV0Helium3PerEvent"), countHe3);
    registry.fill(HIST("hV0AntiHelium3PerEvent"), countAntiHe3);

    for ( auto& v0 : V0s ){
      registry.fill(HIST("hV0McCheckCounter"), 0.5);
      registry.fill(HIST("hV0CollisionID"), v0.collisionId());
      auto recopostrack = v0.posTrack_as<MyTracks>();
      auto reconegtrack = v0.negTrack_as<MyTracks>();
      if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle()) {
        continue;
      }
      registry.fill(HIST("hV0McCheckCounter"), 1.5);
      auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();
      auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
      registry.fill(HIST("hV0McPostrackPt"), mcpostrack.pt());
      registry.fill(HIST("hV0McNegtrackPt"), mcnegtrack.pt());

      registry.fill(HIST("hV0PostrackPt"), v0.posTrack_as<MyTracks>().pt());
      registry.fill(HIST("hV0NegtrackPt"), v0.negTrack_as<MyTracks>().pt());
      registry.fill(HIST("hV0PostrackTPCNClus"), v0.posTrack_as<MyTracks>().tpcNClsCrossedRows());
      registry.fill(HIST("hV0NegtrackTPCNClus"), v0.negTrack_as<MyTracks>().tpcNClsCrossedRows());
      registry.fill(HIST("hV0PostrackDcatoPv"), v0.posTrack_as<MyTracks>().dcaXY());
      registry.fill(HIST("hV0NegtrackDcatoPv"), v0.negTrack_as<MyTracks>().dcaXY());

      registry.fill(HIST("hV0DausTPCBB"), v0.posTrack_as<MyTracks>().p()*v0.posTrack_as<MyTracks>().sign(), v0.posTrack_as<MyTracks>().tpcSignal());
      registry.fill(HIST("hV0DausTPCBB"), v0.negTrack_as<MyTracks>().p()*v0.negTrack_as<MyTracks>().sign(), v0.negTrack_as<MyTracks>().tpcSignal());

      uint32_t pTrackPID = v0.posTrack_as<MyTracks>().pidForTracking();
      uint32_t nTrackPID = v0.negTrack_as<MyTracks>().pidForTracking();
      registry.fill(HIST("hV0PostrackID"), pTrackPID);
      registry.fill(HIST("hV0NegtrackID"), nTrackPID);
      int pTrackCharge = o2::track::pid_constants::sCharges[pTrackPID];
      int nTrackCharge = o2::track::pid_constants::sCharges[nTrackPID];
      registry.fill(HIST("hV0PostrackCharge"), pTrackCharge);
      registry.fill(HIST("hV0NegtrackCharge"), nTrackCharge);

      registry.fill(HIST("hNSigmaHelium3Pos"), GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ) );
      registry.fill(HIST("hNSigmaHelium3Neg"), GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) );
      registry.fill(HIST("hNSigmaPionPos"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaPionNeg"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaTritonPos"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaTritonNeg"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());

      v0mcinfo( mcpostrack.pdgCode(), pTrackPID, v0.posTrack_as<MyTracks>().px(), v0.posTrack_as<MyTracks>().py(), v0.posTrack_as<MyTracks>().pz(), mcpostrack.px(), mcpostrack.py(), mcpostrack.pz(), v0.posTrack_as<MyTracks>().tpcNClsCrossedRows(), v0.posTrack_as<MyTracks>().tpcSignal(), GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ), v0.posTrack_as<MyTracks>().tpcNSigmaPi(), v0.posTrack_as<MyTracks>().tpcNSigmaTr(), v0.posTrack_as<MyTracks>().dcaXY(),  
          mcnegtrack.pdgCode(), nTrackPID, v0.negTrack_as<MyTracks>().px(), v0.negTrack_as<MyTracks>().py(), v0.negTrack_as<MyTracks>().pz(), mcnegtrack.px(), mcnegtrack.py(), mcnegtrack.pz(), v0.negTrack_as<MyTracks>().tpcNClsCrossedRows(), v0.negTrack_as<MyTracks>().tpcSignal(), GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ), v0.negTrack_as<MyTracks>().tpcNSigmaPi(), v0.negTrack_as<MyTracks>().tpcNSigmaTr(), v0.negTrack_as<MyTracks>().dcaXY() );


      if (mcpostrack.pdgCode() == 1000020030 && mcnegtrack.pdgCode() == -211){
        for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
          for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
            if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
              registry.fill(HIST("hTrueHypertritonCounter"), 0.5);
              if (particleMotherOfNeg.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueHypertritonCounter"), 1.5);
              }
            }
          }
        }
      }
      if (mcnegtrack.pdgCode() == -1000020030 && mcpostrack.pdgCode() == 211) {
        for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
          for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
            if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
              registry.fill(HIST("hTrueHypertritonCounter"), 6.5);
              if (particleMotherOfNeg.isPhysicalPrimary()) {
                registry.fill(HIST("hTrueHypertritonCounter"), 7.5);
              }
            }
          }
        }
      }

      auto pTrack = getTrackParCov(v0.posTrack_as<MyTracks>());
      auto nTrack = getTrackParCov(v0.negTrack_as<MyTracks>());
      auto pTrackCopy = o2::track::TrackParCov(pTrack);
      auto nTrackCopy = o2::track::TrackParCov(nTrack);

      std::array<float, 3> pos = {0.};
      std::array<float, 3> pvec0 = {0.};
      std::array<float, 3> pvec1 = {0.};

      int nCand = fitter.process(pTrackCopy, nTrackCopy);
      if (nCand == 0) {
        continue;
      }
      double finalXpos = fitter.getTrack(0).getX();
      double finalXneg = fitter.getTrack(1).getX();

      // Rotate to desired alpha
      pTrack.rotateParam(fitter.getTrack(0).getAlpha());
      nTrack.rotateParam(fitter.getTrack(1).getAlpha());
      o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
      o2::base::Propagator::Instance()->propagateToX(pTrack, finalXpos, d_bz, maxSnp, maxStep, matCorr);
      o2::base::Propagator::Instance()->propagateToX(nTrack, finalXneg, d_bz, maxSnp, maxStep, matCorr);

      nCand = fitter.process(pTrack, nTrack);
      if (nCand == 0) {
        continue;
      }

      pTrack.getPxPyPzGlo(pvec0);
      nTrack.getPxPyPzGlo(pvec1);
      const auto& vtx = fitter.getPCACandidate();

      double hypertritonMcMass = -999, hypertritonMass = -999;
      double hypertritonRapidity = -999, hypertritonPt = -999, V0radius = -999, V0CosinePA = -999, ct = -999;
      float He3Nsigma = -999, PionNSigma = -999, He3Pt = -999, PionPt = -999, PionDcatoPV = -999;
      if (mcpostrack.pdgCode() == 1000020030 && mcnegtrack.pdgCode() == -211){
        for (int i = 0; i < 3; i++) {
          pos[i] = vtx[i];
          pvec0[i] = 2*pvec0[i];
        }
        hypertritonMcMass = RecoDecay::m(array{array{mcpostrack.px(), mcpostrack.py(), mcpostrack.pz()}, array{mcnegtrack.px(), mcnegtrack.py(), mcnegtrack.pz()}}, array{2.80839, RecoDecay::getMassPDG(kPiPlus)}); 
        hypertritonMass = RecoDecay::m(array{array{pvec0[0], pvec0[1], pvec0[2]}, array{pvec1[0], pvec1[1], pvec1[2]}}, array{2.80839, RecoDecay::getMassPDG(kPiPlus)}); 
        He3Nsigma   = GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ) ;
        PionNSigma  = v0.negTrack_as<MyTracks>().tpcNSigmaPi();
        He3Pt       = TMath::Sqrt(pvec0[0]*pvec0[0]+pvec0[1]*pvec0[1]);
        PionPt      = TMath::Sqrt(pvec1[0]*pvec1[0]+pvec1[1]*pvec1[1]);
        PionDcatoPV = v0.negTrack_as<MyTracks>().dcaXY();

        registry.fill(HIST("hV0DauHelium3McPx"), mcpostrack.px());
        registry.fill(HIST("hV0DauHelium3McPy"), mcpostrack.py());
        registry.fill(HIST("hV0DauHelium3McPz"), mcpostrack.pz());
        registry.fill(HIST("hV0DauHelium3McPt"), mcpostrack.pt());
        registry.fill(HIST("hV0DauHelium3McP"), mcpostrack.p());

        registry.fill(HIST("hV0DauHelium3Px"), pvec0[0]);
        registry.fill(HIST("hV0DauHelium3Py"), pvec0[1]);
        registry.fill(HIST("hV0DauHelium3Pz"), pvec0[2]);
        registry.fill(HIST("hV0DauHelium3Pt"), TMath::Sqrt(pvec0[0]*pvec0[0]+pvec0[1]*pvec0[1]) );
        registry.fill(HIST("hV0DauHelium3P"), TMath::Sqrt(pvec0[0]*pvec0[0]+pvec0[1]*pvec0[1]+pvec0[2]*pvec0[2]) );

        registry.fill(HIST("hV0DauHelium3DiffPx"), pvec0[0]-mcpostrack.px());
        registry.fill(HIST("hV0DauHelium3DiffPy"), pvec0[1]-mcpostrack.py());
        registry.fill(HIST("hV0DauHelium3DiffPz"), pvec0[2]-mcpostrack.pz());
        registry.fill(HIST("hV0DauHelium3DiffPt"), TMath::Sqrt(pvec0[0]*pvec0[0]+pvec0[1]*pvec0[1]) - mcpostrack.pt());
        registry.fill(HIST("hV0DauHelium3DiffP"), TMath::Sqrt(pvec0[0]*pvec0[0]+pvec0[1]*pvec0[1]+pvec0[2]*pvec0[2]) - mcpostrack.p());
        registry.fill(HIST("hV0DauHelium3TPCBB"), v0.posTrack_as<MyTracks>().p()*v0.posTrack_as<MyTracks>().sign(), v0.posTrack_as<MyTracks>().tpcSignal());
        registry.fill(HIST("hV0DauHelium3NSigmaHelium3"), GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() ) );
        registry.fill(HIST("hV0DauHelium3NSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
        registry.fill(HIST("hV0DauHelium3NSigmaTriton"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
        registry.fill(HIST("hV0DauPionNSigmaHelium3"), GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) );
        registry.fill(HIST("hV0DauPionNSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
        registry.fill(HIST("hV0McHypertritonMass"), hypertritonMcMass);
        registry.fill(HIST("hV0HypertritonMass"), hypertritonMass);

      }

      if (mcnegtrack.pdgCode() == -1000020030 && mcpostrack.pdgCode() == 211) {
        for (int i = 0; i < 3; i++) {
          pos[i] = vtx[i];
          pvec1[i] = 2*pvec1[i];
        }

        hypertritonMcMass = RecoDecay::m(array{array{mcpostrack.px(), mcpostrack.py(), mcpostrack.pz()}, array{mcnegtrack.px(), mcnegtrack.py(), mcnegtrack.pz()}}, array{RecoDecay::getMassPDG(kPiPlus), 2.80839}); 
        hypertritonMass = RecoDecay::m(array{array{pvec0[0], pvec0[1], pvec0[2]}, array{pvec1[0], pvec1[1], pvec1[2]}}, array{RecoDecay::getMassPDG(kPiPlus), 2.80839}); 
        He3Nsigma   = GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) ;
        PionNSigma  = v0.posTrack_as<MyTracks>().tpcNSigmaPi();
        He3Pt       = TMath::Sqrt(pvec1[0]*pvec1[0]+pvec1[1]*pvec1[1]);
        PionPt      = TMath::Sqrt(pvec0[0]*pvec0[0]+pvec0[1]*pvec0[1]);
        PionDcatoPV = v0.posTrack_as<MyTracks>().dcaXY();


        registry.fill(HIST("hV0DauHelium3McPx"), mcnegtrack.px());
        registry.fill(HIST("hV0DauHelium3McPy"), mcnegtrack.py());
        registry.fill(HIST("hV0DauHelium3McPz"), mcnegtrack.pz());
        registry.fill(HIST("hV0DauHelium3McPt"), mcnegtrack.pt());
        registry.fill(HIST("hV0DauHelium3McP"), mcnegtrack.p());

        registry.fill(HIST("hV0DauHelium3Px"), pvec1[0]);
        registry.fill(HIST("hV0DauHelium3Py"), pvec1[1]);
        registry.fill(HIST("hV0DauHelium3Pz"), pvec1[2]);
        registry.fill(HIST("hV0DauHelium3Pt"), TMath::Sqrt(pvec1[0]*pvec1[0]+pvec1[1]*pvec1[1]) );
        registry.fill(HIST("hV0DauHelium3P"), TMath::Sqrt(pvec1[0]*pvec1[0]+pvec1[1]*pvec1[1]+pvec1[2]*pvec1[2]) );

        registry.fill(HIST("hV0DauHelium3DiffPx"), pvec1[0]-mcpostrack.px());
        registry.fill(HIST("hV0DauHelium3DiffPy"), pvec1[1]-mcpostrack.py());
        registry.fill(HIST("hV0DauHelium3DiffPz"), pvec1[2]-mcpostrack.pz());
        registry.fill(HIST("hV0DauHelium3DiffPt"), TMath::Sqrt(pvec1[0]*pvec1[0]+pvec1[1]*pvec1[1]) - mcpostrack.pt());
        registry.fill(HIST("hV0DauHelium3DiffP"), TMath::Sqrt(pvec1[0]*pvec1[0]+pvec1[1]*pvec1[1]+pvec1[2]*pvec1[2]) - mcpostrack.p());

        registry.fill(HIST("hV0DauHelium3TPCBB"), v0.negTrack_as<MyTracks>().p()*v0.posTrack_as<MyTracks>().sign(), v0.posTrack_as<MyTracks>().tpcSignal());
        registry.fill(HIST("hV0DauHelium3NSigmaHelium3"), GetTPCNSigmaHe3(2*v0.negTrack_as<MyTracks>().p(), v0.negTrack_as<MyTracks>().tpcSignal() ) );
        registry.fill(HIST("hV0DauHelium3NSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
        registry.fill(HIST("hV0DauHelium3NSigmaTriton"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());
        registry.fill(HIST("hV0DauPionNSigmaHelium3"), GetTPCNSigmaHe3(2*v0.posTrack_as<MyTracks>().p(), v0.posTrack_as<MyTracks>().tpcSignal() )  );
        registry.fill(HIST("hV0DauPionNSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
        registry.fill(HIST("hV0McHypertritonMass"), hypertritonMcMass);
        registry.fill(HIST("hV0HypertritonMass"), hypertritonMass);

      }

      hypertritonRapidity = RecoDecay::y(array{pvec0[0]+pvec1[0], pvec0[1]+pvec1[1], pvec0[2]+pvec1[2]}, 2.991); 
      hypertritonPt = TMath::Sqrt(TMath::Power(pvec0[0]+pvec1[0], 2) + TMath::Power(pvec0[1]+pvec1[1], 2));
      V0radius = RecoDecay::sqrtSumOfSquares(pos[0], pos[1]);
      V0CosinePA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
      ct = std::sqrt(std::pow(collision.posX() - pos[0], 2) + std::pow(collision.posY() - pos[1], 2) + std::pow(collision.posZ() - pos[2], 2)) / (TMath::Sqrt( TMath::Power(pvec0[0]+pvec1[0], 2) + TMath::Power( pvec0[1]+pvec1[1], 2) + TMath::Power( pvec0[2]+pvec1[2], 2)) + 1E-10) * 2.991;

      if (v0.posTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows && v0.negTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows){
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass1"), hypertritonMass);
      if (fitter.getChi2AtPCACandidate() > dcav0dau) {
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass2"), hypertritonMass);
      if (V0CosinePA < v0cospa) {
        continue;
      }
      registry.fill(HIST("hV0McHypertritonMassAfterV0Cut"), hypertritonMcMass);
      registry.fill(HIST("hV0HypertritonMass3"), hypertritonMass);
      if (V0radius > v0radius) {
        registry.fill(HIST("hV0HypertritonMassTest"), hypertritonMass);
      }
      for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
        for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
          if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
            registry.fill(HIST("hTrueHypertritonCounter"), 2.5);
            if (particleMotherOfNeg.isPhysicalPrimary()) {
              registry.fill(HIST("hTrueHypertritonCounter"), 3.5);
            }
          }
        }
      }
      for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
        for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
          if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
            registry.fill(HIST("hTrueHypertritonCounter"), 8.5);
            if (particleMotherOfNeg.isPhysicalPrimary()) {
              registry.fill(HIST("hTrueHypertritonCounter"), 9.5);
            }
          }
        }
      }

      if (TMath::Abs(He3Nsigma) > TpcPidNsigmaCut){
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass4"), hypertritonMass);
      if (TMath::Abs(PionNSigma) > TpcPidNsigmaCut ) {
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass5"), hypertritonMass);
      if ( TMath::Abs(v0.posTrack_as<MyTracks>().eta()) > etacut || TMath::Abs(v0.negTrack_as<MyTracks>().eta()) > etacut){
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass6"), hypertritonMass);
      if (TMath::Abs(hypertritonRapidity) > rapidity){ 
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass7"), hypertritonMass);
      if (ct > lifetimecut) {
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass8"), hypertritonMass);
      if ( PionPt < 0.2 || PionPt > 1.2 ) {
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass9"), hypertritonMass);
      if ( He3Pt < 1.8 || He3Pt > 10 ) {
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass10"), hypertritonMass);
      if (hypertritonPt < 2 || hypertritonPt > 9 ) {
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass11"), hypertritonMass);
      if (TMath::Abs(PionDcatoPV) < dcapiontopv) {
        continue;
      }
      registry.fill(HIST("hV0HypertritonMass12"), hypertritonMass);

      for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
        for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
          if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
            registry.fill(HIST("hTrueHypertritonCounter"), 4.5);
            if (particleMotherOfNeg.isPhysicalPrimary()) {
              registry.fill(HIST("hTrueHypertritonCounter"), 5.5);
            }
          }
        }
      }
      for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
        for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
          if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
            registry.fill(HIST("hTrueHypertritonCounter"), 10.5);
            if (particleMotherOfNeg.isPhysicalPrimary()) {
              registry.fill(HIST("hTrueHypertritonCounter"), 11.5);
            }
          }
        }
      }

    }
  }
};


struct V0DataInitializer {
  Spawns<aod::V0Datas> v0datas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertritonAnalysisMc>(cfgc),
      adaptAnalysisTask<hypertritonQa>(cfgc),
      adaptAnalysisTask<hypertritonMcParticleCount>(cfgc),
      adaptAnalysisTask<hypertritonTrackMcinfo>(cfgc),
      adaptAnalysisTask<goodtrackcheck>(cfgc),
      adaptAnalysisTask<V0DataInitializer>(cfgc),
      adaptAnalysisTask<V0McCheck>(cfgc)
  };
}
