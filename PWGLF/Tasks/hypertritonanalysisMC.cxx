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
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    aimeric.landou@cern.ch (MC adaptation)
//    david.dobrigkeit.chinellato@cern.ch (original hypertritonanalysis task)

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

#include "DetectorsVertexing/DCAFitterN.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>
#include "DetectorsBase/Propagator.h"

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

        {"hSelectedV0Counter", "hSelectedV0Counter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},

        {"hHypertritonFeedDownMatrix", "hHypertritonFeedDownMatrix", {HistType::kTH2F, {{200, 0.0f, 10.0f, "#it{p}_{T}^{#Hypertriton} (GeV/c)"}, {200, 0.0f, 10.0f, "#it{p}_{T}^{#Omega-} (GeV/c)"}}}},
        {"hAntiHypertritonFeedDownMatrix", "hAntiHypertritonFeedDownMatrix", {HistType::kTH2F, {{200, 0.0f, 10.0f, "#it{p}_{T}^{#bar{#Hypertriton}} (GeV/c)"}, {200, 0.0f, 10.0f, "#it{p}_{T}^{#Omega+} (GeV/c)"}}}},

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
  Configurable<float> etacut{"etacut", 0.9, "etacut"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcapiontopv{"dcapiontopv", .1, "DCA Pion To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.9, "rapidity"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  //Configurable<bool> boolArmenterosCut{"boolArmenterosCut", true, "cut on Armenteros-Podolanski graph"};
  //Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  //Configurable<bool> hasItsTest{"hasItsTest", false, "hasItsTest"};

  static constexpr float defaultLifetimeCuts[1][1] = {{40.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 1, {"lifetimecutHypertriton"}}, "lifetimecut"};

  Filter preFilterV0 = aod::v0data::dcaV0daughters < dcav0dau;

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
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
      if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {
        registry.fill(HIST("hSelectedV0Counter"), 7.5);

        if(v0.negTrack_as<MyTracks>().pt() > 0.2 && v0.negTrack_as<MyTracks>().pt() < 1.2 && v0.posTrack_as<MyTracks>().pt() > 1.8 && v0.posTrack_as<MyTracks>().pt() < 10 && v0.pt() > 2 && v0.pt() < 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 8.5);
          if (TMath::Abs(v0.dcanegtopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 9.5);

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
      if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut ) {

        registry.fill(HIST("hSelectedV0Counter"), 7.5);
        if(v0.posTrack_as<MyTracks>().pt() > 0.2 && v0.posTrack_as<MyTracks>().pt() < 1.2 && v0.negTrack_as<MyTracks>().pt() > 1.8 && v0.negTrack_as<MyTracks>().pt() < 10 && v0.pt() < 2 && v0.pt() > 9 ){
          registry.fill(HIST("hSelectedV0Counter"), 8.5);
          if (TMath::Abs(v0.dcapostopv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedV0Counter"), 9.5);

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
      registry.fill(HIST("hSelectedV0Counter"), 1.5);
      if(v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("hSelectedV0Counter"), 2.5);

        auto reconegtrack = v0.negTrack_as<MyTracks>();
        auto recopostrack = v0.posTrack_as<MyTracks>();

        auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
        auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();

        if (TMath::Abs(v0.yHypertriton()) < rapidity) {
          registry.fill(HIST("hSelectedV0Counter"), 3.5);
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991 < lifetimecut->get("lifetimecutHypertriton")) {
            registry.fill(HIST("hSelectedV0Counter"), 4.5);
            registry.fill(HIST("h3dMassHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mHypertriton());
            registry.fill(HIST("h3dMassAntiHypertriton"), collision.centRun2V0M(), v0.pt(), v0.mAntiHypertriton());
            registry.fill(HIST("hArmenterosPostAnalyserCuts"), v0.alpha(), v0.qtarm());

            for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
              for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
                if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
                  registry.fill(HIST("h3dMassHypertriton_MC_truePt"), collision.centRun2V0M(), particleMotherOfNeg.pt(), v0.mHypertriton());
                  registry.fill(HIST("hArmenterosPostAnalyserCuts_MC"), v0.alpha(), v0.qtarm());
                }
                if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -1010010030) {
                  registry.fill(HIST("h3dMassAntiHypertriton_MC_truePt"), collision.centRun2V0M(), particleMotherOfNeg.pt(), v0.mAntiHypertriton());
                  registry.fill(HIST("hArmenterosPostAnalyserCuts_MC"), v0.alpha(), v0.qtarm());
                }

                if (saveDcaHist == 1) {
                  registry.fill(HIST("h3dMassHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mHypertriton());
                  registry.fill(HIST("h3dMassAntiHypertritonDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiHypertriton());

                  if (particleMotherOfNeg.isPhysicalPrimary() && particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 1010010030) {
                    registry.fill(HIST("h3dMassHypertritonDca_MC_truePt"), v0.dcaV0daughters(), particleMotherOfNeg.pt(), v0.mHypertriton());
                  }
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


struct hypertritonTrackCount {
  //Basic checks
  HistogramRegistry registry{
    "registry",
      {

        {"hTotalCollCounter", "hTotalCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hParticleCount", "hParticleCount", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hParticleCount2", "hParticleCount2", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hTPCNClsCrossedRows", "hTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
        {"hTestCount", "hHelium3CountBeforeCut", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hTrackEta", "hTrackEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hTrackMcRapidity", "hTrackMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hTrackNsigmaHelium3", "hTrackNsigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hTrackNsigmaPion", "hTrackNsigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3Eta", "hHelium3Eta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHelium3McRapidity", "hHelium3McRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHypertritonEta", "hHypertritomEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHypertritonMcRapidity", "hHypertritonMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},


        {"hHypertritonMcPt", "hHypertritonMcPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},

        {"hHelium3Px", "hHelium3Px", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3Py", "hHelium3Py", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3Pz", "hHelium3Pz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3Pt", "hHelium3Pt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hHelium3P", "hHelium3P", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hHelium3McPx", "hHelium3McPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3McPy", "hHelium3McPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3McPz", "hHelium3McPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hHelium3McPt", "hHelium3McPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hHelium3McP", "hHelium3McP", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hPionPx", "hPionPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionPy", "hPionPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionPz", "hPionPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionPt", "hPionPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hPionP", "hPionP", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hPionMcPx", "hPionMcPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionMcPy", "hPionMcPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionMcPz", "hPionMcPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hPionMcPt", "hPionMcPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hPionMcP", "hPionMcP", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},


        {"hHelium3NsigmaHelium3", "hHelium3NsigmaHelium3", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3NsigmaPion", "hHelium3NsigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3NsigmaTriton", "hHelium3NsigmaTriton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3TPCBB", "hHelium3TPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal" }}}},
        {"hHelium3NsigmaHelium3Test", "hHelium3NsigmaHelium3Test", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
        {"hHelium3TPCBBTest", "hHelium3TPCBBTest", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal" }}}},
        {"hHelium3WrongNSigmaDCAToPV", "hHelium3WrongNSigmaDCAToPV", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
        {"hHelium3WrongNSigmaPt", "hHelium3WrongNSigmaPt", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},
        {"hHelium3WrongNSigmaP", "hHelium3WrongNSigmaP", {HistType::kTH1F, {{150, 0.0f, 15.0f}}}},

        {"hHelium3TPCNClsCrossedRows", "hHelium3TPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
        {"hPionTPCNClsCrossedRows", "hPionTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},

        {"hTPCBB", "hTPCBB", {HistType::kTH2F, {{120, -8.0f, 8.0f, "p/z(GeV/c)"}, {100, 0.0f, 1000.0f, "TPCSignal" }}}},
      },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(2, "Has_mcparticle");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(3, "Rapidity Cut(off)");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(4, "McisHelium3");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(5, "McisHypertriton");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(6, "McisPion");

    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(2, "Has_mcparticle");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(3, "Rapidity Cut(off)");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(4, "McisHelium3");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(5, "McisHypertriton");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(6, "McisPion");
  }


  void process(aod::Collision const& collision, aod::McParticles const& mcparticles, MyTracks const& tracks)
  {

    registry.fill(HIST("hTotalCollCounter"), 0.5);

    for (auto& track : tracks) {
      registry.fill(HIST("hParticleCount"), 0.5);
      if (track.tpcNClsCrossedRows() > 70){
        registry.fill(HIST("hParticleCount2"), 0.5);
      }
      registry.fill(HIST("hTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("hTrackNsigmaHelium3"), track.tpcNSigmaHe());
      registry.fill(HIST("hTrackNsigmaPion"), track.tpcNSigmaPi());

      if(!track.has_mcParticle()){
        continue;
      }
      registry.fill(HIST("hParticleCount"), 1.5);
      if (track.tpcNClsCrossedRows() > 70){
        registry.fill(HIST("hParticleCount2"), 1.5);
      }
      auto mcparticle = track.mcParticle_as<aod::McParticles>();
      registry.fill(HIST("hTPCBB"), track.p()*track.sign(), track.tpcSignal());
      if (mcparticle.pdgCode() == 1000020030 || mcparticle.pdgCode() == -1000020030) {
        registry.fill(HIST("hTestCount"), 0.5);
      }

      //if (TMath::Abs(mcparticle.y()) > 0.9) {continue;}
      registry.fill(HIST("hParticleCount"), 2.5);
      if (track.tpcNClsCrossedRows() > 70){
        registry.fill(HIST("hParticleCount2"), 2.5);
      }
      registry.fill(HIST("hTrackEta"), track.eta());
      registry.fill(HIST("hTrackMcRapidity"), mcparticle.y());


      if (mcparticle.pdgCode() == 1000020030 || mcparticle.pdgCode() == -1000020030) {
        registry.fill(HIST("hParticleCount"), 3.5);
        if (track.tpcNClsCrossedRows() > 70){
          registry.fill(HIST("hParticleCount2"), 3.5);
          registry.fill(HIST("hHelium3NsigmaHelium3Test"), track.tpcNSigmaHe());
          registry.fill(HIST("hHelium3TPCBBTest"), track.p()*track.sign(), track.tpcSignal());
          if (TMath::Abs(track.tpcNSigmaHe()) > 5){
            registry.fill(HIST("hHelium3WrongNSigmaDCAToPV"), track.dcaXY());
            registry.fill(HIST("hHelium3WrongNSigmaPt"), 2*track.pt());
            registry.fill(HIST("hHelium3WrongNSigmaP"), 2*track.p());
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

        registry.fill(HIST("hHelium3NsigmaHelium3"), track.tpcNSigmaHe());
        registry.fill(HIST("hHelium3NsigmaPion"), track.tpcNSigmaPi());
        registry.fill(HIST("hHelium3NsigmaTriton"), track.tpcNSigmaTr());
        registry.fill(HIST("hHelium3TPCNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hHelium3Eta"), track.eta());
        registry.fill(HIST("hHelium3McRapidity"), mcparticle.y());
        registry.fill(HIST("hHelium3TPCBB"), track.p()*track.sign(), track.tpcSignal());

      }

      if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hParticleCount"), 4.5);
        if (track.tpcNClsCrossedRows() > 70){
          registry.fill(HIST("hParticleCount2"), 4.5);
        }
        registry.fill(HIST("hHypertritonMcPt"), mcparticle.pt());
        registry.fill(HIST("hHypertritonEta"), track.eta());
        registry.fill(HIST("hHypertritonMcRapidity"), mcparticle.y());
      }

      if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
        registry.fill(HIST("hParticleCount"), 5.5);
        if (track.tpcNClsCrossedRows() > 70){
          registry.fill(HIST("hParticleCount2"), 5.5);
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


struct hypertritonParticleCountMc {
  //Basic checks
  Produces<aod::McHypertritonCheck> hypertrtitonCheckTable;
  Produces<aod::McHelium3Check> helium3CheckTable;
  HistogramRegistry registry{
    "registry",
      {
        {"hTotalMcCollCounter", "hTotalMcCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hSelAndRecoMcCollCounter", "hSelAndRecoMcCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},

        {"hMcHypertritonEffCheck", "hMcHypertritonCount", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
        {"hMcParticleCount", "hMcParticleCount", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hHypertritonCount_PtDiff", "hHypertritonCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hAntiHypertritonCount_PtDiff", "hAntiHypertritonCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},


        {"hMcHypertritonCount", "hMcHypertritonCount", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
        {"hMcAntiHypertritonCount", "hMcAntiHypertritonCount", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},

        {"hMcHypertritonCheck", "hMcHypertritonCheck", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
        {"hMcHelium3Check", "hMcHelium3Check", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},

        {"hMcHypertritonPt", "hMcHypertritonPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hMcHelium3Pt", "hMcHelium3Pt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hMcPionPt", "hMcPionPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hMcHypertritonEffCheck"))->GetXaxis()->SetBinLabel(1, "PdgCode");
    registry.get<TH1>(HIST("hMcHypertritonEffCheck"))->GetXaxis()->SetBinLabel(2, "PtCut");
    registry.get<TH1>(HIST("hMcHypertritonEffCheck"))->GetXaxis()->SetBinLabel(3, "Lifetime(off)");
    registry.get<TH1>(HIST("hMcHypertritonEffCheck"))->GetXaxis()->SetBinLabel(4, "Rapidity");

    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(1, "ReadinAfterEvtRecSelCut(off)");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(2, "y<0.9(off)");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(3, "IsPhysicalPrimary");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(4, "(Anti)Helium3");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(5, "(Anti)Hypertriton");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(6, "HasDaughter");
    registry.get<TH1>(HIST("hMcHypertritonCount"))->GetXaxis()->SetBinLabel(1, "primary Hypertriton mothers");
    registry.get<TH1>(HIST("hMcHypertritonCount"))->GetXaxis()->SetBinLabel(2, "decaying into V0");
    registry.get<TH1>(HIST("hMcAntiHypertritonCount"))->GetXaxis()->SetBinLabel(1, "primary AntiHypertriton mothers");
    registry.get<TH1>(HIST("hMcAntiHypertritonCount"))->GetXaxis()->SetBinLabel(2, "decaying into V0");
    registry.get<TH1>(HIST("hMcHelium3Check"))->GetXaxis()->SetBinLabel(1, "PdgCode");
    registry.get<TH1>(HIST("hMcHelium3Check"))->GetXaxis()->SetBinLabel(2, "has_mothers");
    registry.get<TH1>(HIST("hMcHelium3Check"))->GetXaxis()->SetBinLabel(3, "has_daughters");
    registry.get<TH1>(HIST("hMcHelium3Check"))->GetXaxis()->SetBinLabel(4, "has_mothers&&has_daughters");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(1, "PdgCode");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(2, "has_mothers");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(3, "has_daughters");
  }

  Configurable<float> rapidityMCcut{"rapidityMCcut", 0.9, "rapidity cut MC count"};
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
        if (mcparticle.has_mothers()){
          registry.fill(HIST("hMcHelium3Check"), 1.5);
        }
        if (mcparticle.has_daughters()){
          registry.fill(HIST("hMcHelium3Check"), 2.5);
        }
        if (mcparticle.has_mothers() && mcparticle.has_daughters()){
          registry.fill(HIST("hMcHelium3Check"), 3.5);
        }
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
        registry.fill(HIST("hMcHypertritonCheck"), 0.5);
        registry.fill(HIST("hMcHypertritonEffCheck"), 0.5);
        registry.fill(HIST("hMcHypertritonPt"), mcparticle.pt());
        if (mcparticle.has_mothers()){
          registry.fill(HIST("hMcHypertritonCheck"), 1.5);
        }
        if (mcparticle.has_daughters()){
          registry.fill(HIST("hMcHypertritonCheck"), 2.5);
        }
        for (auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          hypertrtitonCheckTable(mcparticleDaughter.pdgCode());
        }
        if (mcparticle.pt() > 2 && mcparticle.pt() < 10){
          registry.fill(HIST("hMcHypertritonEffCheck"), 1.5);
          //if (mcparticle.Lifetime() < 40){
          registry.fill(HIST("hMcHypertritonEffCheck"), 2.5);
          if (TMath::Abs(mcparticle.y()) < 0.8){
            registry.fill(HIST("hMcHypertritonEffCheck"), 3.5);
          }
          //}
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

      if (mcparticle.pdgCode() == 1010010030) {
        registry.fill(HIST("hMcHypertritonCount"), 0.5);
        for (auto& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
          for (auto& mcparticleDaughter1 : mcparticle.daughters_as<aod::McParticles>()) {
            if (mcparticleDaughter0.pdgCode() == -211 && mcparticleDaughter1.pdgCode() == 1000020030) {
              registry.fill(HIST("hMcHypertritonCount"), 1.5);
              registry.fill(HIST("hHypertritonCount_PtDiff"), mcparticle.pt());
            }
          }
        }
      }
      if (mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hMcAntiHypertritonCount"), 0.5);
        for (auto& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
          for (auto& mcparticleDaughter1 : mcparticle.daughters_as<aod::McParticles>()) {
            if (mcparticleDaughter0.pdgCode() == 211 && mcparticleDaughter1.pdgCode() == -1000020030) {
              registry.fill(HIST("hMcAntiHypertritonCount"), 1.5);
              registry.fill(HIST("hAntiHypertritonCount_PtDiff"), mcparticle.pt());
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
    DECLARE_SOA_COLUMN(PosTrackPx, posTrackPx, double);
    DECLARE_SOA_COLUMN(PosTrackPy, posTrackPy, double);
    DECLARE_SOA_COLUMN(PosTrackPz, posTrackPz, double);
    DECLARE_SOA_COLUMN(PosTrackMcPx, posTrackMcPx, double);
    DECLARE_SOA_COLUMN(PosTrackMcPy, posTrackMcPy, double);
    DECLARE_SOA_COLUMN(PosTrackMcPz, posTrackMcPz, double);
    DECLARE_SOA_COLUMN(PosTrackDiffPx, posTrackDiffPx, double);
    DECLARE_SOA_COLUMN(PosTrackDiffPy, posTrackDiffPy, double);
    DECLARE_SOA_COLUMN(PosTrackDiffPz, posTrackDiffPz, double);
    DECLARE_SOA_COLUMN(PosTrackTPCNCls, posTrackTPCNCls, int);
    DECLARE_SOA_COLUMN(PosTrackTPCSignals, posTrackTPCSignals, double);
    DECLARE_SOA_COLUMN(PosTrackTPCNSigmaHelium3, posTrackTPCNSigmaHe, double);
    DECLARE_SOA_COLUMN(PosTrackTPCNSigmaPion, posTrackTPCNSigmaPi, double);
    DECLARE_SOA_COLUMN(PosTrackTPCNSigmaTriton, posTrackTPCNSigmaTr, double);
    DECLARE_SOA_COLUMN(PosTrackDcaXY, posTrackDcaXY, double);

    DECLARE_SOA_COLUMN(NegTrackMcPdgCode, negTrackMcPdgCode, int);
    DECLARE_SOA_COLUMN(NegTrackPdgPID, negTrackPdgPID, int);
    DECLARE_SOA_COLUMN(NegTrackPx, negTrackPx, double);
    DECLARE_SOA_COLUMN(NegTrackPy, negTrackPy, double);
    DECLARE_SOA_COLUMN(NegTrackPz, negTrackPz, double);
    DECLARE_SOA_COLUMN(NegTrackMcPx, negTrackMcPx, double);
    DECLARE_SOA_COLUMN(NegTrackMcPy, negTrackMcPy, double);
    DECLARE_SOA_COLUMN(NegTrackMcPz, negTrackMcPz, double);
    DECLARE_SOA_COLUMN(NegTrackDiffPx, negTrackDiffPx, double);
    DECLARE_SOA_COLUMN(NegTrackDiffPy, negTrackDiffPy, double);
    DECLARE_SOA_COLUMN(NegTrackDiffPz, negTrackDiffPz, double);
    DECLARE_SOA_COLUMN(NegTrackTPCNCls, negTrackTPCNCls, int);
    DECLARE_SOA_COLUMN(NegTrackTPCSignals, negTrackTPCSignals, int);
    DECLARE_SOA_COLUMN(NegTrackTPCNSigmaHelium3, negTrackTPCNSigmaHe, double);
    DECLARE_SOA_COLUMN(NegTrackTPCNSigmaPion, negTrackTPCNSigmaPion, double);
    DECLARE_SOA_COLUMN(NegTrackTPCNSigmaTriton, negTrackTPCNSigmaTriton, double);
    DECLARE_SOA_COLUMN(NegTrackDcaXY, negTrackDcaXY, double);
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
      v0mc::PosTrackDiffPx,
      v0mc::PosTrackDiffPy,
      v0mc::PosTrackDiffPz,
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
      v0mc::NegTrackDiffPx, 
      v0mc::NegTrackDiffPy, 
      v0mc::NegTrackDiffPz, 
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
        {"hV0PostrackPt", "hV0PostrackPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hV0NegtrackPt", "hV0PostrackPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hV0PostrackTPCNClus", "hV0PostrackTPCNClus", {HistType::kTH1F, {{24, 0.0f, 240.0f}}}},
        {"hV0NegtrackTPCNClus", "hV0NegtrackTPCNClus", {HistType::kTH1F, {{24, 0.0f, 240.0f}}}},
        {"hV0PostrackDcatoPv", "hV0PostrackDcatoPv", {HistType::kTH1F, {{10, 0.0f, 1.0f}}}},
        {"hV0NegtrackDcatoPv", "hV0NegtrackDcatoPv", {HistType::kTH1F, {{10, 0.0f, 1.0f}}}},
        {"hV0McPostrackPt", "hV0McPostrackPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hV0McNegtrackPt", "hV0McNegtrackPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
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
        {"hV0DauHelium3Pt", "hV0DauHelium3Pt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hV0DauHelium3P", "hV0DauHelium3P", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hV0DauHelium3McPx", "hV0DauHelium3McPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3McPy", "hV0DauHelium3McPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3McPz", "hV0DauHelium3McPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3McPt", "hV0DauHelium3McPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hV0DauHelium3McP", "hV0DauHelium3McP", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPx", "hV0DauHelium3DiffPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPy", "hV0DauHelium3DiffPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPz", "hV0DauHelium3DiffPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffPt", "hV0DauHelium3DiffPt", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3DiffP", "hV0DauHelium3DiffP", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
        {"hV0DauHelium3TPCBB", "hV0DauHelium3TPCBB", {HistType::kTH2F, {{120, -8.0f, 8.0f, "p/z(GeV/c)"}, {100, 0.0f, 1000.0f, "TPCSignal" }}}},
        /*{"hV0PionPx", "hV0PionPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionPy", "hV0PionPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionPz", "hV0PionPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionPt", "hV0PionPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
          {"hV0PionP", "hV0PionP", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
          {"hV0PionMcPx", "hV0PionMcPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionMcPy", "hV0PionMcPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionMcPz", "hV0PionMcPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionMcPt", "hV0PionMcPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
          {"hV0PionMcP", "hV0PionMcP", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
          {"hV0PionDiffPx", "hV0PionDiffPx", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionDiffPy", "hV0PionDiffPy", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionDiffPz", "hV0PionDiffPz", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionDiffPt", "hV0PionDiffPt", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},
          {"hV0PionDiffP", "hV0PionDiffP", {HistType::kTH1F, {{400, -10.0f, 10.0f}}}},*/
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

      registry.fill(HIST("hNSigmaHelium3Pos"), v0.posTrack_as<MyTracks>().tpcNSigmaHe());
      registry.fill(HIST("hNSigmaHelium3Neg"), v0.negTrack_as<MyTracks>().tpcNSigmaHe());
      registry.fill(HIST("hNSigmaPionPos"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaPionNeg"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
      registry.fill(HIST("hNSigmaTritonPos"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
      registry.fill(HIST("hNSigmaTritonNeg"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());

      v0mcinfo( mcpostrack.pdgCode(), pTrackPID, v0.posTrack_as<MyTracks>().px(), v0.posTrack_as<MyTracks>().py(), v0.posTrack_as<MyTracks>().pz(), mcpostrack.px(), mcpostrack.py(), mcpostrack.pz(), v0.posTrack_as<MyTracks>().px() - mcpostrack.px(), v0.posTrack_as<MyTracks>().py() - mcpostrack.py(), v0.posTrack_as<MyTracks>().pz() - mcpostrack.pz(), v0.posTrack_as<MyTracks>().tpcNClsCrossedRows(), v0.posTrack_as<MyTracks>().tpcSignal(),v0.posTrack_as<MyTracks>().tpcNSigmaHe(), v0.posTrack_as<MyTracks>().tpcNSigmaPi(), v0.posTrack_as<MyTracks>().tpcNSigmaTr(), v0.posTrack_as<MyTracks>().dcaXY(),  
          mcnegtrack.pdgCode(), nTrackPID, v0.negTrack_as<MyTracks>().px(), v0.negTrack_as<MyTracks>().py(), v0.negTrack_as<MyTracks>().pz(), mcnegtrack.px(), mcnegtrack.py(), mcnegtrack.pz(), v0.negTrack_as<MyTracks>().px() - mcnegtrack.px(), v0.negTrack_as<MyTracks>().py() - mcnegtrack.py(), v0.negTrack_as<MyTracks>().pz() - mcnegtrack.pz(),  v0.negTrack_as<MyTracks>().tpcNClsCrossedRows(),  v0.negTrack_as<MyTracks>().tpcSignal(), v0.negTrack_as<MyTracks>().tpcNSigmaHe(), v0.negTrack_as<MyTracks>().tpcNSigmaPi(), v0.negTrack_as<MyTracks>().tpcNSigmaTr(), v0.negTrack_as<MyTracks>().dcaXY() );


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

        bool fillV0Hist = true;
        auto pTrack = getTrackParCov(v0.posTrack_as<MyTracks>());
        auto nTrack = getTrackParCov(v0.negTrack_as<MyTracks>());
        auto pTrackCopy = o2::track::TrackParCov(pTrack);
        auto nTrackCopy = o2::track::TrackParCov(nTrack);

        std::array<float, 3> pos = {0.};
        std::array<float, 3> pvec0 = {0.};
        std::array<float, 3> pvec1 = {0.};

        int nCand = fitter.process(pTrackCopy, nTrackCopy);
        if (nCand == 0) {
          fillV0Hist = false;
        }
        if(fillV0Hist){
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
            fillV0Hist = false;
          }

          pTrack.getPxPyPzGlo(pvec0);
          nTrack.getPxPyPzGlo(pvec1);
          const auto& vtx = fitter.getPCACandidate();
          for (int i = 0; i < 3; i++) {
            pos[i] = vtx[i];
            pvec0[i] = 2*pvec0[i];
          }
        }
        double hypertritonMcMass = RecoDecay::m(array{array{mcpostrack.px(), mcpostrack.py(), mcpostrack.pz()}, array{mcnegtrack.px(), mcnegtrack.py(), mcnegtrack.pz()}}, array{2.80839, RecoDecay::getMassPDG(kPiPlus)}); 
        double hypertritonMass = RecoDecay::m(array{array{pvec0[0], pvec0[1], pvec0[2]}, array{pvec1[0], pvec1[1], pvec1[2]}}, array{2.80839, RecoDecay::getMassPDG(kPiPlus)}); 
        double hypertritonRapidity = RecoDecay::y(array{pvec0[0]+pvec1[0], pvec0[1]+pvec1[1], pvec0[2]+pvec1[2]}, 2.991); 
        double hypertritonPt = TMath::Sqrt(TMath::Power(pvec0[0]+pvec1[1], 2) + TMath::Power(pvec0[1]+pvec1[1], 2));
        auto V0radius = RecoDecay::sqrtSumOfSquares(pos[0], pos[1]);
        auto V0CosinePA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
        double ct = std::sqrt(std::pow(collision.posX() - pos[0], 2) + std::pow(collision.posY() - pos[1], 2) + std::pow(collision.posZ() - pos[2], 2)) / (TMath::Sqrt( TMath::Power(pvec0[0]+pvec1[0], 2) + TMath::Power( pvec0[1]+pvec1[1], 2) + TMath::Power( pvec0[2]+pvec1[2], 2)) + 1E-10) * 2.991;

        if (fillV0Hist){

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
          registry.fill(HIST("hV0DauHelium3NSigmaHelium3"), v0.posTrack_as<MyTracks>().tpcNSigmaHe());
          registry.fill(HIST("hV0DauHelium3NSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
          registry.fill(HIST("hV0DauHelium3NSigmaTriton"), v0.posTrack_as<MyTracks>().tpcNSigmaTr());
          registry.fill(HIST("hV0DauPionNSigmaHelium3"), v0.negTrack_as<MyTracks>().tpcNSigmaHe());
          registry.fill(HIST("hV0DauPionNSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
          registry.fill(HIST("hV0McHypertritonMass"), hypertritonMcMass);
          registry.fill(HIST("hV0HypertritonMass"), hypertritonMass);
        }



        for (int i=0; i<1; i++){
          if (v0.posTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows && v0.negTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows){
            break;
          }
          registry.fill(HIST("hV0HypertritonMass1"), hypertritonMass);
          if (fitter.getChi2AtPCACandidate() > dcav0dau) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass2"), hypertritonMass);
          if (V0CosinePA < v0cospa) {
            break;
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
          if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaHe()) > TpcPidNsigmaCut){
            break;
          }
          registry.fill(HIST("hV0HypertritonMass4"), hypertritonMass);
          if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) > TpcPidNsigmaCut ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass5"), hypertritonMass);
          if ( TMath::Abs(v0.posTrack_as<MyTracks>().eta()) > etacut || TMath::Abs(v0.negTrack_as<MyTracks>().eta()) > etacut){
            continue;
          }
          registry.fill(HIST("hV0HypertritonMass6"), hypertritonMass);
          if (hypertritonRapidity > rapidity){ 
            break;
          }
          registry.fill(HIST("hV0HypertritonMass7"), hypertritonMass);
          if (ct > lifetimecut) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass8"), hypertritonMass);
          if (v0.negTrack_as<MyTracks>().pt() < 0.2 || v0.negTrack_as<MyTracks>().pt() > 1.2 ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass9"), hypertritonMass);
          if (v0.posTrack_as<MyTracks>().pt() < 1.8 || v0.posTrack_as<MyTracks>().pt() > 10 ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass10"), hypertritonMass);
          if (hypertritonPt < 2 || hypertritonPt > 9 ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass11"), hypertritonMass);
          if (TMath::Abs(v0.negTrack_as<MyTracks>().dcaXY()) < dcapiontopv) {
            break;
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

        bool fillV0Hist = true;
        auto pTrack = getTrackParCov(v0.posTrack_as<MyTracks>());
        auto nTrack = getTrackParCov(v0.negTrack_as<MyTracks>());
        auto pTrackCopy = o2::track::TrackParCov(pTrack);
        auto nTrackCopy = o2::track::TrackParCov(nTrack);

        std::array<float, 3> pos = {0.};
        std::array<float, 3> pvec0 = {0.};
        std::array<float, 3> pvec1 = {0.};

        int nCand = fitter.process(pTrackCopy, nTrackCopy);
        if (nCand == 0) {
          fillV0Hist = false;
        }
        if(fillV0Hist){
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
            fillV0Hist = false;
          }

          pTrack.getPxPyPzGlo(pvec0);
          nTrack.getPxPyPzGlo(pvec1);
          const auto& vtx = fitter.getPCACandidate();
          for (int i = 0; i < 3; i++) {
            pos[i] = vtx[i];
            pvec1[i] = 2*pvec1[i];
          }
        }

        double hypertritonMcMass = RecoDecay::m(array{array{mcpostrack.px(), mcpostrack.py(), mcpostrack.pz()}, array{mcnegtrack.px(), mcnegtrack.py(), mcnegtrack.pz()}}, array{RecoDecay::getMassPDG(kPiPlus), 2.80839}); 
        double hypertritonMass = RecoDecay::m(array{array{pvec0[0], pvec0[1], pvec0[2]}, array{pvec1[0], pvec1[1], pvec1[2]}}, array{RecoDecay::getMassPDG(kPiPlus), 2.80839}); 
        double hypertritonRapidity = RecoDecay::y(array{pvec0[0]+pvec1[0], pvec0[1]+pvec1[1], pvec0[2]+pvec1[2]}, 2.991); 
        double hypertritonPt = TMath::Sqrt(TMath::Power(pvec0[0]+pvec1[1], 2) + TMath::Power(pvec0[1]+pvec1[1], 2));
        auto V0radius = RecoDecay::sqrtSumOfSquares(pos[0], pos[1]);
        auto V0CosinePA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
        double ct = std::sqrt(std::pow(collision.posX() - pos[0], 2) + std::pow(collision.posY() - pos[1], 2) + std::pow(collision.posZ() - pos[2], 2)) / (TMath::Sqrt( TMath::Power(pvec0[0]+pvec1[0], 2) + TMath::Power( pvec0[1]+pvec1[1], 2) + TMath::Power( pvec0[2]+pvec1[2], 2)) + 1E-10) * 2.991;

        if (fillV0Hist){

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
          registry.fill(HIST("hV0DauHelium3NSigmaHelium3"), v0.negTrack_as<MyTracks>().tpcNSigmaHe());
          registry.fill(HIST("hV0DauHelium3NSigmaPion"), v0.negTrack_as<MyTracks>().tpcNSigmaPi());
          registry.fill(HIST("hV0DauHelium3NSigmaTriton"), v0.negTrack_as<MyTracks>().tpcNSigmaTr());
          registry.fill(HIST("hV0DauPionNSigmaHelium3"), v0.posTrack_as<MyTracks>().tpcNSigmaHe());
          registry.fill(HIST("hV0DauPionNSigmaPion"), v0.posTrack_as<MyTracks>().tpcNSigmaPi());
          registry.fill(HIST("hV0McHypertritonMass"), hypertritonMcMass);
          registry.fill(HIST("hV0HypertritonMass"), hypertritonMass);
        }



        for (int i=0; i<1; i++){
          if (v0.posTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows && v0.negTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows){
            break;
          }
          registry.fill(HIST("hV0HypertritonMass1"), hypertritonMass);
          if (fitter.getChi2AtPCACandidate() > dcav0dau) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass2"), hypertritonMass);
          if (V0CosinePA < v0cospa) {
            break;
          }
          registry.fill(HIST("hV0McHypertritonMassAfterV0Cut"), hypertritonMcMass);
          registry.fill(HIST("hV0HypertritonMass3"), hypertritonMass);
          if (V0radius > v0radius) {
            registry.fill(HIST("hV0HypertritonMassTest"), hypertritonMass);
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
          if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaHe()) > TpcPidNsigmaCut){
            break;
          }
          registry.fill(HIST("hV0HypertritonMass4"), hypertritonMass);
          if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) > TpcPidNsigmaCut ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass5"), hypertritonMass);
          if ( TMath::Abs(v0.posTrack_as<MyTracks>().eta()) > etacut || TMath::Abs(v0.negTrack_as<MyTracks>().eta()) > etacut){
            continue;
          }
          registry.fill(HIST("hV0HypertritonMass6"), hypertritonMass);
          if (hypertritonRapidity > rapidity){ 
            break;
          }
          registry.fill(HIST("hV0HypertritonMass7"), hypertritonMass);
          if (ct > lifetimecut) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass8"), hypertritonMass);
          if (v0.posTrack_as<MyTracks>().pt() < 0.2 || v0.posTrack_as<MyTracks>().pt() > 1.2 ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass9"), hypertritonMass);
          if (v0.negTrack_as<MyTracks>().pt() < 1.8 || v0.negTrack_as<MyTracks>().pt() > 10 ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass10"), hypertritonMass);
          if (hypertritonPt < 2 || hypertritonPt > 9 ) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass11"), hypertritonMass);
          if (TMath::Abs(v0.posTrack_as<MyTracks>().dcaXY()) < dcapiontopv) {
            break;
          }
          registry.fill(HIST("hV0HypertritonMass12"), hypertritonMass);

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
      adaptAnalysisTask<hypertritonParticleCountMc>(cfgc),
      adaptAnalysisTask<hypertritonTrackCount>(cfgc),
      adaptAnalysisTask<V0DataInitializer>(cfgc),
      adaptAnalysisTask<V0McCheck>(cfgc)};
}
