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
#include "../Hypv0Table.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using MyTracks = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullHe, aod::pidTPCFullTr, aod::pidTPCFullKa, aod::pidTPCFullPr>;


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
      },
  };

  void init(InitContext const&)
  {
    AxisSpec massAxisHypertriton = {300, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {massAxisHypertriton}});
    registry.add("hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {massAxisHypertriton}});
  }

  void process(aod::Collision const& collision, aod::HypV0Datas const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
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
        {"MCmomID_Hypertriton", "MCmomID_Hypertriton", {HistType::kTH1I, {{4000000, 0, 4000000}}}},
        {"MCmomID_AntiHypertriton", "MCmomID_AntiHypertriton", {HistType::kTH1I, {{4000000, 0, 4000000}}}},
        {"V0loopFiltersCounts", "V0loopFiltersCounts", {HistType::kTH1F, {{11, 0.0f, 11.0f}}}},

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

  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisHypertriton = {massHypertritonbinning, "Inv. Mass (GeV/c^{2})"};

    registry.add("h3dMassHypertritonDca", "h3dMassHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
    registry.add("h3dMassAntiHypertritonDca", "h3dMassAntiHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
    registry.add("h3dMassHypertritonDca_MC_truePt", "h3dMassHypertritonDca_MC_truePt", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
    registry.add("h3dMassAntiHypertritonDca_MC_truePt", "h3dMassAntiHypertritonDca_MC_truePt", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});

    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(1, "V0 Candidates");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(2, "V0Radius and CosPA");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(4, "Hypertriton Rapidity");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(5, "Hypertriton lifetime cut");
    registry.get<TH1>(HIST("V0loopFiltersCounts"))->GetXaxis()->SetBinLabel(6, "Hypertriton TPC PID cut");
  }

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcapiontopv{"dcapiontopv", .1, "DCA Pion To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.9, "rapidity"};
  Configurable<int> saveDcaHist{"saveDcaHist", 0, "saveDcaHist"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 3, "TpcPidNsigmaCut"};
  //Configurable<bool> boolArmenterosCut{"boolArmenterosCut", true, "cut on Armenteros-Podolanski graph"};
  //Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  //Configurable<bool> hasItsTest{"hasItsTest", false, "hasItsTest"};

  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutHypertriton", "lifetimecutK0S"}}, "lifetimecut"};

  Filter preFilterV0 = aod::hypv0data::dcaV0daughters < dcav0dau;

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::HypV0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
    // void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
  {
    /*if (eventSelection && !collision.sel8()) {
      return;
      }*/
    registry.fill(HIST("hSelectedEventCounter"), 0.5);

    for (auto& v0 : fullV0s) {
      //   FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      registry.fill(HIST("V0loopFiltersCounts"), 0.5);
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("V0loopFiltersCounts"), 1.5);

        auto reconegtrack = v0.negTrack_as<MyTracks>();
        auto recopostrack = v0.posTrack_as<MyTracks>();
        if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle()) {
          continue;
        }
        auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
        auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();

        if (TMath::Abs(v0.yHypertriton()) < rapidity) {
          registry.fill(HIST("V0loopFiltersCounts"), 3.5);
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991 < lifetimecut->get("lifetimecutHypertriton")) {
            registry.fill(HIST("V0loopFiltersCounts"), 4.5);

            //Hypertriton
            if (TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(v0.dcanegtopv()) > dcapiontopv) { //previous 900Gev pp analysis had nSigma< 5 for pt<0.7Gev and tpcNSigmaStorePr<3 for pt>0.7GeV; and no cut on K0S
              registry.fill(HIST("V0loopFiltersCounts"), 5.5);
              // registry.fill(HIST("h3dMassHypertriton"), collision.centV0M(), v0.pt(), v0.mHypertriton());
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

            // AntiHypertriton
            if (TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaHe()) < TpcPidNsigmaCut && TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(v0.dcapostopv()) > dcapiontopv) { //previous 900Gev pp analysis had nSigma< 5 for pt<0.7Gev and tpcNSigmaStorePr<3 for pt>0.7GeV; and no cut on K0S
              registry.fill(HIST("V0loopFiltersCounts"), 5.5);
              // registry.fill(HIST("h3dMassHypertriton"), collision.centV0M(), v0.pt(), v0.mHypertriton());
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
  }
  PROCESS_SWITCH(hypertritonAnalysisMc, processRun3, "Process Run 3 data", true);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<aod::HypV0Datas> const& fullV0s, aod::McParticles const& mcParticles, MyTracks const& tracks)
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
      registry.fill(HIST("V0loopFiltersCounts"), 0.5);
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("V0loopFiltersCounts"), 1.5);

        auto reconegtrack = v0.negTrack_as<MyTracks>();
        auto recopostrack = v0.posTrack_as<MyTracks>();

        auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
        auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();

        if (TMath::Abs(v0.yHypertriton()) < rapidity) {
          registry.fill(HIST("V0loopFiltersCounts"), 3.5);
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * 2.991 < lifetimecut->get("lifetimecutHypertriton")) {
            registry.fill(HIST("V0loopFiltersCounts"), 4.5);
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

struct hypertritonParticleCountMc {
  //Basic checks
  HistogramRegistry registry{
    "registry",
      {
        {"hMcParticleCount", "hMcParticleCount", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
        {"hHypertritonCount", "hHypertritonCount", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
        {"hAntiHypertritonCount", "hAntiHypertritonCount", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
        {"hHypertritonCount_PtDiff", "hHypertritonCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
        {"hAntiHypertritonCount_PtDiff", "hAntiHypertritonCount_PtDiff", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},

        {"hSelAndRecoMcCollCounter", "hSelAndRecoMcCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hTotalMcCollCounter", "hTotalMcCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},

      },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(2, "y<0.9");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(3, "IsPhysicalPrimary");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(4, "(Anti)Helium3");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(5, "(Anti)Hypertriton");
    registry.get<TH1>(HIST("hMcParticleCount"))->GetXaxis()->SetBinLabel(6, "HasDaughter");
    registry.get<TH1>(HIST("hHypertritonCount"))->GetXaxis()->SetBinLabel(1, "primary Hypertriton mothers");
    registry.get<TH1>(HIST("hHypertritonCount"))->GetXaxis()->SetBinLabel(2, "decaying into V0");
    registry.get<TH1>(HIST("hAntiHypertritonCount"))->GetXaxis()->SetBinLabel(1, "primary AntiHypertriton mothers");
    registry.get<TH1>(HIST("hAntiHypertritonCount"))->GetXaxis()->SetBinLabel(2, "decaying into V0");
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

    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();

    registry.fill(HIST("hTotalMcCollCounter"), 0.5);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    registry.fill(HIST("hSelAndRecoMcCollCounter"), 0.5);

    for (auto& mcparticle : mcParticles) {
      registry.fill(HIST("hMcParticleCount"), 0.5);
      if (TMath::Abs(mcparticle.y()) < rapidityMCcut) {
        registry.fill(HIST("hMcParticleCount"), 1.5);
        if (mcparticle.isPhysicalPrimary()) {
          registry.fill(HIST("hMcParticleCount"), 2.5);
        }
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
            registry.fill(HIST("hHypertritonCount"), 0.5);
            for (auto& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
              for (auto& mcparticleDaughter1 : mcparticle.daughters_as<aod::McParticles>()) {
                if (mcparticleDaughter0.pdgCode() == -211 && mcparticleDaughter1.pdgCode() == 10000200030) {
                  registry.fill(HIST("hHypertritonCount"), 1.5);
                  registry.fill(HIST("hHypertritonCount_PtDiff"), mcparticle.pt());
                }
              }
            }
          }
          if (mcparticle.pdgCode() == -1010010030) {
            registry.fill(HIST("hAntiHypertritonCount"), 0.5);
            for (auto& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
              for (auto& mcparticleDaughter1 : mcparticle.daughters_as<aod::McParticles>()) {
                if (mcparticleDaughter0.pdgCode() == 211 && mcparticleDaughter1.pdgCode() == -1000020030) {
                  registry.fill(HIST("hAntiHypertritonCount"), 1.5);
                  registry.fill(HIST("hAntiHypertritonCount_PtDiff"), mcparticle.pt());
                }
              }
            }
          }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertritonAnalysisMc>(cfgc),
      adaptAnalysisTask<hypertritonQa>(cfgc),
      adaptAnalysisTask<hypertritonParticleCountMc>(cfgc)};
}
