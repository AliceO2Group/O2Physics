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
// Example StoredVtx3BodyDatas analysis task
// ========================
//
// This code loops over a StoredVtx3BodyDatas table and produces some
// standard analysis output. It requires either
// the hypertriton3bodybuilder or hypertriton3bodyfinder (not recommaended) tasks
// to have been executed in the workflow (before).
//
// author: yuanzhe.wang@cern.ch

#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCDe, aod::pidTPCPr>;
using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::pidTPCFullPr, aod::McTrackLabels>;

struct hypertriton3bodyQaMc {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxRadius", "hVtxRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hVtxCosPA", "hVtxCosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hDCATrack0ToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCATrack1ToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, 10.0f, 10.0f, "cm"}}}},
      {"hDCATrack2ToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, 10.0f, 10.0f, "cm"}}}},
      {"hDCAVtxDau", "hDCAVtxDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hVtxPt", "hVtxPt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTrack0Pt", "hTrack0Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTrack1Pt", "hTrack1Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTrack2Pt", "hTrack2Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
    },
  };
  void init(InitContext const&)
  {
    AxisSpec massAxis = {120, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {massAxis}});
    registry.add("hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {massAxis}});
  }
  void process(aod::Collision const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, MyTracks const& tracks)
  {

    for (auto& vtx : vtx3bodydatas) {
      registry.fill(HIST("hVtxRadius"), vtx.vtxradius());
      registry.fill(HIST("hVtxCosPA"), vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCATrack0ToPV"), vtx.dcatrack0topv());
      registry.fill(HIST("hDCATrack1ToPV"), vtx.dcatrack1topv());
      registry.fill(HIST("hDCATrack2ToPV"), vtx.dcatrack2topv());
      registry.fill(HIST("hDCAVtxDau"), vtx.dcaVtxdaughters());
      registry.fill(HIST("hVtxPt"), vtx.pt());
      registry.fill(HIST("hTrack0Pt"), vtx.track0pt());
      registry.fill(HIST("hTrack1Pt"), vtx.track1pt());
      registry.fill(HIST("hTrack2Pt"), vtx.track2pt());
      registry.fill(HIST("hMassHypertriton"), vtx.mHypertriton());
      registry.fill(HIST("hMassAntiHypertriton"), vtx.mAntiHypertriton());
    }
  }
};

struct hypertriton3bodyAnalysisMc {

  HistogramRegistry registry{
    "registry",
    {
      {"hSelectedEventCounter", "hSelectedEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hSelectedCandidatesCounter", "hSelectedCandidatesCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
      {"hSelectedTrueHypertritonCounter", "hSelectedTrueHypertritonCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
      {"hTestCounter", "hTestCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
      {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{100, 2.95f, 3.05f}}}},
      {"hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {{100, 2.95f, 3.05f}}}},
      {"hMassHypertritonTotal", "hMassHypertritonTotal", {HistType::kTH1F, {{300, 2.9f, 3.2f}}}},
      {"hPtProton", "hPtProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiPion", "hPtAntiPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtDeuteron", "hPtDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiProton", "hPtAntiProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtPion", "hPtPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiDeuteron", "hPtAntiDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"h3dMassHypertriton", "h3dMassHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dMassAntiHypertriton", "h3dMassAntiHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dTotalHypertriton", "h3dTotalHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dTotalTrueHypertriton", "h3dTotalTrueHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
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

    if (saveDcaHist == 1) {
      registry.add("h3dMassHypertritonDca", "h3dMassHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
      registry.add("h3dMassAntiHypertritonDca", "h3dMassAntiHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
    }

    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(2, "VtxCosPA");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(3, "TrackEta");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(4, "MomRapidity");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(5, "Lifetime");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(6, "DcaV0Dau");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(7, "TPCPID");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(8, "PtCut");
    registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(9, "PionDcatoPV");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(2, "VtxCosPA");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(3, "TrackEta");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(4, "MomRapidity");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(5, "Lifetime");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(6, "DcaV0Dau");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(7, "TPCPID");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(8, "PtCut");
    registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(9, "PionDcatoPV");
  }

  // Selection criteria
  Configurable<double> vtxcospa{"vtxcospa", 0.9, "Vtx CosPA"};          // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"}; // loose cut
  Configurable<float> dcapiontopv{"dcapiontopv", .00, "DCA Pion To PV"};
  Configurable<float> etacut{"etacut", 1, "etacut"};
  Configurable<float> rapidity{"rapidity", 1, "rapidity"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<float> lifetimecut{"lifetimecut", 40., "lifetimecut"};
  Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
  Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
  Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
  Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
  Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
  Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, MyTracks const& tracks, aod::McParticles const& particlesMC)
  {
    registry.fill(HIST("hSelectedEventCounter"), 0.5);
    if (eventSelection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hSelectedEventCounter"), 1.5);

    for (auto& vtx : vtx3bodydatas) {
      // couldn't filter cosPA and radius variables (dynamic columns)

      // int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      double MClifetime = -1;
      bool isFromHypertriton = false;
      auto lTrack0 = vtx.track0_as<MyTracks>();
      auto lTrack1 = vtx.track1_as<MyTracks>();
      auto lTrack2 = vtx.track2_as<MyTracks>();
      if (lTrack0.has_mcParticle() && lTrack1.has_mcParticle() && lTrack2.has_mcParticle()) {
        auto lMCTrack0 = lTrack0.mcParticle_as<aod::McParticles>();
        auto lMCTrack1 = lTrack1.mcParticle_as<aod::McParticles>();
        auto lMCTrack2 = lTrack2.mcParticle_as<aod::McParticles>();
        if (lMCTrack0.has_mothers() && lMCTrack1.has_mothers() && lMCTrack2.has_mothers()) {
          for (auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
            for (auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
              for (auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
                if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
                  // lLabel = lMother1.globalIndex();
                  lPt = lMother1.pt();
                  lPDG = lMother1.pdgCode();
                  if ((lPDG == 1010010030 && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == 1000010020) ||
                      (lPDG == -1010010030 && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -1000010020)) {
                    isFromHypertriton = true;
                    MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
                  }
                }
              }
            }
          }
        }
      }

      registry.fill(HIST("hSelectedCandidatesCounter"), 0.5);
      if (isFromHypertriton) {
        registry.fill(HIST("hSelectedTrueHypertritonCounter"), 0.5);
      }

      if (vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()) < vtxcospa) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), 1.5);
      if (isFromHypertriton) {
        registry.fill(HIST("hSelectedTrueHypertritonCounter"), 1.5);
      }
      if (TMath::Abs(vtx.track0_as<MyTracks>().eta()) > etacut || TMath::Abs(vtx.track1_as<MyTracks>().eta()) > etacut || TMath::Abs(vtx.track2_as<MyTracks>().eta()) > etacut) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), 2.5);
      if (isFromHypertriton) {
        registry.fill(HIST("hSelectedTrueHypertritonCounter"), 2.5);
      }
      if (TMath::Abs(vtx.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), 3.5);
      if (isFromHypertriton) {
        registry.fill(HIST("hSelectedTrueHypertritonCounter"), 3.5);
      }
      double ct = vtx.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassHyperTriton;
      if (ct > lifetimecut) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), 4.5);
      if (isFromHypertriton) {
        registry.fill(HIST("hSelectedTrueHypertritonCounter"), 4.5);
      }
      if (vtx.dcaVtxdaughters() > dcavtxdau) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), 5.5);
      if (isFromHypertriton) {
        registry.fill(HIST("hSelectedTrueHypertritonCounter"), 5.5);
      }

      // 3sigma region for Dalitz plot
      double mcsigma = 0.0015; // from local MC
      double lowerlimit = o2::constants::physics::MassHyperTriton - 3 * mcsigma;
      double upperlimit = o2::constants::physics::MassHyperTriton + 3 * mcsigma;

      // Hypertriton
      if (TMath::Abs(vtx.track0_as<MyTracks>().tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(vtx.track1_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(vtx.track2_as<MyTracks>().tpcNSigmaDe()) < TpcPidNsigmaCut) {
        registry.fill(HIST("hSelectedCandidatesCounter"), 6.5);
        if (isFromHypertriton) {
          registry.fill(HIST("hSelectedTrueHypertritonCounter"), 6.5);
        }

        registry.fill(HIST("hTestCounter"), 0.5);
        if (vtx.track1pt() > minPionPt && vtx.track1pt() < maxPionPt) {
          registry.fill(HIST("hTestCounter"), 1.5);
          if (vtx.track0pt() > minProtonPt && vtx.track0pt() < maxProtonPt) {
            registry.fill(HIST("hTestCounter"), 2.5);
            if (vtx.pt() > minDeuteronPt && vtx.pt() < maxDeuteronPt) {
              registry.fill(HIST("hTestCounter"), 3.5);
            }
          }
        }

        if (vtx.track0pt() > minProtonPt && vtx.track0pt() < maxProtonPt && vtx.track1pt() > minPionPt && vtx.track1pt() < maxPionPt && vtx.track2pt() > minDeuteronPt && vtx.track2pt() < maxDeuteronPt) {
          registry.fill(HIST("hSelectedCandidatesCounter"), 7.5);
          if (isFromHypertriton) {
            registry.fill(HIST("hSelectedTrueHypertritonCounter"), 7.5);
          }

          if (TMath::Abs(vtx.dcatrack1topv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedCandidatesCounter"), 8.5);
            if (isFromHypertriton) {
              registry.fill(HIST("hSelectedTrueHypertritonCounter"), 8.5);
            }

            registry.fill(HIST("hPtProton"), vtx.track0pt());
            registry.fill(HIST("hPtAntiPion"), vtx.track1pt());
            registry.fill(HIST("hPtDeuteron"), vtx.track2pt());
            registry.fill(HIST("hMassHypertriton"), vtx.mHypertriton());
            registry.fill(HIST("hMassHypertritonTotal"), vtx.mHypertriton());
            registry.fill(HIST("h3dMassHypertriton"), 0., vtx.pt(), vtx.mHypertriton()); // collision.centV0M() instead of 0. once available
            registry.fill(HIST("h3dTotalHypertriton"), ct, vtx.pt(), vtx.mHypertriton());
            if (vtx.mHypertriton() > lowerlimit && vtx.mHypertriton() < upperlimit) {
              registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{vtx.pxtrack0(), vtx.pytrack0(), vtx.pztrack0()}, array{vtx.pxtrack2(), vtx.pytrack2(), vtx.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}),
                            RecoDecay::m2(array{array{vtx.pxtrack0(), vtx.pytrack0(), vtx.pztrack0()}, array{vtx.pxtrack1(), vtx.pytrack1(), vtx.pztrack1()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
            }

            if (isFromHypertriton) {
              registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, vtx.mHypertriton());
            }
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassHypertritonDca"), vtx.dcaVtxdaughters(), vtx.pt(), vtx.mHypertriton());
            }
          }
        }
      }

      // AntiHypertriton
      if (TMath::Abs(vtx.track0_as<MyTracks>().tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(vtx.track1_as<MyTracks>().tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(vtx.track2_as<MyTracks>().tpcNSigmaDe()) < TpcPidNsigmaCut) {

        registry.fill(HIST("hSelectedCandidatesCounter"), 6.5);
        if (isFromHypertriton) {
          registry.fill(HIST("hSelectedTrueHypertritonCounter"), 6.5);
        }

        registry.fill(HIST("hTestCounter"), 0.5);
        if (vtx.track0pt() > minPionPt && vtx.track0pt() < maxPionPt) {
          registry.fill(HIST("hTestCounter"), 1.5);
          if (vtx.track1pt() > minProtonPt && vtx.track1pt() < maxProtonPt) {
            registry.fill(HIST("hTestCounter"), 2.5);
            if (vtx.pt() > minDeuteronPt && vtx.pt() < maxDeuteronPt) {
              registry.fill(HIST("hTestCounter"), 3.5);
            }
          }
        }

        if (vtx.track1pt() > minPionPt && vtx.track1pt() < maxPionPt && vtx.track0pt() > minProtonPt && vtx.track0pt() < maxProtonPt && vtx.track2pt() > minDeuteronPt && vtx.track2pt() < maxDeuteronPt) {
          registry.fill(HIST("hSelectedCandidatesCounter"), 7.5);
          if (isFromHypertriton) {
            registry.fill(HIST("hSelectedTrueHypertritonCounter"), 7.5);
          }
          if (TMath::Abs(vtx.dcatrack0topv()) > dcapiontopv) {
            registry.fill(HIST("hSelectedCandidatesCounter"), 8.5);
            if (isFromHypertriton) {
              registry.fill(HIST("hSelectedTrueHypertritonCounter"), 8.5);
            }

            registry.fill(HIST("hPtAntiProton"), vtx.track1pt());
            registry.fill(HIST("hPtPion"), vtx.track0pt());
            registry.fill(HIST("hPtAntiDeuteron"), vtx.track2pt());
            registry.fill(HIST("hMassAntiHypertriton"), vtx.mAntiHypertriton());
            registry.fill(HIST("hMassHypertritonTotal"), vtx.mAntiHypertriton());
            registry.fill(HIST("h3dMassAntiHypertriton"), 0., vtx.pt(), vtx.mAntiHypertriton()); // collision.centV0M() instead of 0. once available
            registry.fill(HIST("h3dTotalHypertriton"), ct, vtx.pt(), vtx.mAntiHypertriton());
            if (vtx.mAntiHypertriton() > lowerlimit && vtx.mAntiHypertriton() < upperlimit) {
              registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{vtx.pxtrack1(), vtx.pytrack1(), vtx.pztrack1()}, array{vtx.pxtrack2(), vtx.pytrack2(), vtx.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}),
                            RecoDecay::m2(array{array{vtx.pxtrack1(), vtx.pytrack1(), vtx.pztrack1()}, array{vtx.pxtrack0(), vtx.pytrack0(), vtx.pztrack0()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
            }

            if (isFromHypertriton) {
              registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, vtx.mAntiHypertriton());
            }
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassAntiHypertritonDca"), vtx.dcaVtxdaughters(), vtx.pt(), vtx.mAntiHypertriton());
            }
          }
        }
      }
    }
  }
};

// check vtx3body with mclabels
struct hypertriton3bodyLabelCheck {
  HistogramRegistry registry{
    "registry",
    {
      {"hLabeledVtxCounter", "hLabeledVtxCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hMassTrueH3L", "hMassTrueH3L", {HistType::kTH1F, {{100, 2.95f, 3.05f}}}},
      {"hMassTrueH3LMatter", "hMassTrueH3LMatter", {HistType::kTH1F, {{100, 2.95f, 3.05f}}}},
      {"hMassTrueH3LAntiMatter", "hMassTrueH3LAntiMatter", {HistType::kTH1F, {{100, 2.95f, 3.05f}}}},
    },
  };

  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("hLabeledVtxCounter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hLabeledVtxCounter"))->GetXaxis()->SetBinLabel(2, "TrueMCH3L");
  }

  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::Vtx3BodyDatas, aod::McVtx3BodyLabels> const& vtx3bodydatas, MyTracks const& tracks, aod::McParticles const& particlesMC)
  {
    if (eventSelection && !collision.sel8()) {
      return;
    }

    for (auto& vtx : vtx3bodydatas) {
      registry.fill(HIST("hLabeledVtxCounter"), 0.5);
      if (vtx.mcParticleId() != -1) {
        auto mcparticle = vtx.mcParticle_as<aod::McParticles>();
        if (mcparticle.pdgCode() == 1010010030) {
          registry.fill(HIST("hLabeledVtxCounter"), 1.5);
          registry.fill(HIST("hMassTrueH3L"), vtx.mHypertriton());
          registry.fill(HIST("hMassTrueH3LMatter"), vtx.mHypertriton());
        } else if (mcparticle.pdgCode() == -1010010030) {
          registry.fill(HIST("hLabeledVtxCounter"), 1.5);
          registry.fill(HIST("hMassTrueH3L"), vtx.mAntiHypertriton());
          registry.fill(HIST("hMassTrueH3LAntiMatter"), vtx.mAntiHypertriton());
        }
      }
    }
  }
};

// check the properties of daughters candidates and true daughters
struct hypertriton3bodyTrackMcinfo {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hTotalCollCounter", "hTotalCollCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      {"hParticleCount", "hParticleCount", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}},
      {"hParticleCount2", "hParticleCount2", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}}, // for tpcncls > 70
      {"hProtonCount", "hProtonCount", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
      {"hPionCount", "hPionCount", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}},
      {"hDeuteronCount", "hDeuteronCount", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
      {"hBackgroundPionCount", "hBackgroundPionCount", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},

      {"hTPCNClsCrossedRows", "hTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
      {"hTrackEta", "hTrackEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hTrackITSNcls", "hTrackITSNcls", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
      {"hTrackMcRapidity", "hTrackMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hTrackNsigmaProton", "hTrackNsigmaProton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hTrackNsigmaPion", "hTrackNsigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hTrackNsigmaDeuteron", "hTrackNsigmaDeuteron", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},

      {"hHypertritonEta", "hHypertritomEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hHypertritonMcRapidity", "hHypertritonMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hHypertritonMcPt", "hHypertritonMcPt", {HistType::kTH1F, {{300, 0.0f, 15.0f}}}},
      {"hProtonEta", "hProtonEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hProtonMcRapidity", "hProtonMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hPionEta", "hPionEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hPionMcRapidity", "hPionMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hDeuteronEta", "hDeuteronEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hDeuteronMcRapidity", "hDeuteronMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},

      {"hDauProtonPt", "hDauProtonPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauProtonMcPt", "hDauProtonMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauProtonNsigmaProton", "hDauProtonNsigmaProton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hDauProtonNsigmaProtonAfterMomCut", "hDauProtonNsigmaProtonAfterMomCut", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hProtonPt", "hProtonPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonP", "hProtonP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonMcPt", "hProtonMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonMcP", "hProtonMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonNsigmaProton", "hProtonNsigmaProton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hProtonTPCNClsCrossedRows", "hProtonTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
      {"hProtonTPCBB", "hProtonTPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hProtonTPCBBAfterTPCNclsCut", "hProtonTPCBBAfterTPCNclsCut", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},

      {"hDauPionPt", "hDauPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauPionMcPt", "hDauPionMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauPionNsigmaPion", "hDauPionNsigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hDauPionNsigmaPionAfterMomCut", "hDauPionNsigmaPionAfterMomCut", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hPionPt", "hPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionP", "hPionP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionMcPt", "hPionMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionMcP", "hPionMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionNsigmaPion", "hPionNsigmaPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hPionTPCNClsCrossedRows", "hPionTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
      {"hPionTPCBB", "hPionTPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hPionTPCBBAfterTPCNclsCut", "hPionTPCBBAfterTPCNclsCut", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},

      {"hDauDeuteronPt", "hDauDeuteronPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauDeuteronMcPt", "hDauDeuteronMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauDeuteronNsigmaDeuteron", "hDauDeuteronNsigmaDeuteron", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hDauDeuteronNsigmaDeuteronAfterMomCut", "hDauDeuteronNsigmaDeuteronAfterMomCut", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hDeuteronPt", "hDeuteronPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronP", "hDeuteronP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronMcPt", "hDeuteronMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronMcP", "hDeuteronMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronNsigmaDeuteron", "hDeuteronNsigmaDeuteron", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hDeuteronTPCNClsCrossedRows", "hDeuteronTPCNClsCrossedRows", {HistType::kTH1F, {{240, 0.0f, 240.0f}}}},
      {"hDeuteronTPCBB", "hDeuteronTPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDeuteronTPCBBAfterTPCNclsCut", "hDeuteronTPCBBAfterTPCNclsCut", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},

      {"hTPCBB", "hTPCBB", {HistType::kTH2F, {{120, -8.0f, 8.0f, "p/z(GeV/c)"}, {100, 0.0f, 1000.0f, "TPCSignal"}}}},
    },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(2, "Has_mcparticle");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(3, "Rapidity Cut(off)");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(4, "McisHypertriton");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(5, "McisProton");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(6, "McisPion");
    registry.get<TH1>(HIST("hParticleCount"))->GetXaxis()->SetBinLabel(7, "McisDeuteron");

    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(2, "Has_mcparticle");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(3, "Rapidity Cut(off)");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(4, "McisHypertriton");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(5, "McisProton");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(6, "McisPion");
    registry.get<TH1>(HIST("hParticleCount2"))->GetXaxis()->SetBinLabel(7, "McisDeuteron");

    registry.get<TH1>(HIST("hProtonCount"))->GetXaxis()->SetBinLabel(1, "hasMom");
    registry.get<TH1>(HIST("hProtonCount"))->GetXaxis()->SetBinLabel(2, "FromHypertriton");
    registry.get<TH1>(HIST("hProtonCount"))->GetXaxis()->SetBinLabel(3, "TPCNcls");
    registry.get<TH1>(HIST("hProtonCount"))->GetXaxis()->SetBinLabel(4, "Eta");
    registry.get<TH1>(HIST("hProtonCount"))->GetXaxis()->SetBinLabel(5, "Pt");
    registry.get<TH1>(HIST("hProtonCount"))->GetXaxis()->SetBinLabel(6, "TPCPID");
    registry.get<TH1>(HIST("hPionCount"))->GetXaxis()->SetBinLabel(1, "hasMom");
    registry.get<TH1>(HIST("hPionCount"))->GetXaxis()->SetBinLabel(2, "FromHypertriton");
    registry.get<TH1>(HIST("hPionCount"))->GetXaxis()->SetBinLabel(3, "TPCNcls");
    registry.get<TH1>(HIST("hPionCount"))->GetXaxis()->SetBinLabel(4, "Eta");
    registry.get<TH1>(HIST("hPionCount"))->GetXaxis()->SetBinLabel(5, "Pt");
    registry.get<TH1>(HIST("hPionCount"))->GetXaxis()->SetBinLabel(6, "TPCPID");
    registry.get<TH1>(HIST("hPionCount"))->GetXaxis()->SetBinLabel(7, "DcatoPV");
    registry.get<TH1>(HIST("hDeuteronCount"))->GetXaxis()->SetBinLabel(1, "hasMom");
    registry.get<TH1>(HIST("hDeuteronCount"))->GetXaxis()->SetBinLabel(2, "FromHypertriton");
    registry.get<TH1>(HIST("hDeuteronCount"))->GetXaxis()->SetBinLabel(3, "TPCNcls");
    registry.get<TH1>(HIST("hDeuteronCount"))->GetXaxis()->SetBinLabel(4, "Eta");
    registry.get<TH1>(HIST("hDeuteronCount"))->GetXaxis()->SetBinLabel(5, "Pt");
    registry.get<TH1>(HIST("hDeuteronCount"))->GetXaxis()->SetBinLabel(6, "TPCPID");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(1, "NotFromHypertriton");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(2, "TPCNcls");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(3, "Eta");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(4, "Pt");
    registry.get<TH1>(HIST("hBackgroundPionCount"))->GetXaxis()->SetBinLabel(5, "TPCPID");
  }

  Configurable<float> dcapiontopv{"dcapiontopv", .05, "DCA Pion To PV"};
  Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
  Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
  Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
  Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
  Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
  Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};

  void process(aod::Collision const& collision, aod::McParticles const& mcparticles, MyTracks const& tracks)
  {

    registry.fill(HIST("hTotalCollCounter"), 0.5);

    for (auto& track : tracks) {
      registry.fill(HIST("hParticleCount"), 0.5);
      registry.fill(HIST("hTrackITSNcls"), track.itsNCls());
      if (track.tpcNClsCrossedRows() > 70) {
        registry.fill(HIST("hParticleCount2"), 0.5);
      }
      registry.fill(HIST("hTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("hTrackNsigmaDeuteron"), track.tpcNSigmaDe());
      registry.fill(HIST("hTrackNsigmaProton"), track.tpcNSigmaPr());
      registry.fill(HIST("hTrackNsigmaPion"), track.tpcNSigmaPi());

      if (!track.has_mcParticle()) {
        continue;
      }
      registry.fill(HIST("hParticleCount"), 1.5);
      auto mcparticle = track.mcParticle_as<aod::McParticles>();
      if (track.tpcNClsCrossedRows() > 70) {
        registry.fill(HIST("hParticleCount2"), 1.5);
      }
      registry.fill(HIST("hTPCBB"), track.p() * track.sign(), track.tpcSignal());

      // if (TMath::Abs(mcparticle.y()) > 0.9) {continue;}
      registry.fill(HIST("hParticleCount"), 2.5);
      if (track.tpcNClsCrossedRows() > 70) {
        registry.fill(HIST("hParticleCount2"), 2.5);
      }
      registry.fill(HIST("hTrackEta"), track.eta());
      registry.fill(HIST("hTrackMcRapidity"), mcparticle.y());

      // Hypertriton
      if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hParticleCount"), 3.5);
        if (track.tpcNClsCrossedRows() > 70) {
          registry.fill(HIST("hParticleCount2"), 3.5);
        }
        registry.fill(HIST("hHypertritonMcPt"), mcparticle.pt());
        registry.fill(HIST("hHypertritonEta"), track.eta());
        registry.fill(HIST("hHypertritonMcRapidity"), mcparticle.y());
      }

      // Proton
      if (mcparticle.pdgCode() == 2212 || mcparticle.pdgCode() == -2212) {
        registry.fill(HIST("hParticleCount"), 4.5);
        if (track.tpcNClsCrossedRows() > 70) {
          registry.fill(HIST("hParticleCount2"), 4.5);
          registry.fill(HIST("hProtonTPCBBAfterTPCNclsCut"), track.p() * track.sign(), track.tpcSignal());
        }

        if (mcparticle.has_mothers()) {
          registry.fill(HIST("hProtonCount"), 0.5);
          for (auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
            if (particleMother.pdgCode() != 1010010030 && mcparticle.pdgCode() != -1010010030) {
              continue;
            }
            bool haveProton = false, havePion = false, haveDeuteron = false;
            bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
            for (auto& mcparticleDaughter : particleMother.daughters_as<aod::McParticles>()) {
              if (mcparticleDaughter.pdgCode() == 2212) {
                haveProton = true;
              }
              if (mcparticleDaughter.pdgCode() == -2212) {
                haveAntiProton = true;
              }
              if (mcparticleDaughter.pdgCode() == 211) {
                havePion = true;
              }
              if (mcparticleDaughter.pdgCode() == -211) {
                haveAntiPion = true;
              }
              if (mcparticleDaughter.pdgCode() == 1000010020) {
                haveDeuteron = true;
              }
              if (mcparticleDaughter.pdgCode() == -1000010020) {
                haveAntiDeuteron = true;
              }
            }
            if (!(haveProton && haveAntiPion && haveDeuteron && particleMother.pdgCode() == 1010010030) && !(haveAntiProton && havePion && haveAntiDeuteron && particleMother.pdgCode() == -1010010030)) {
              continue;
            }
            registry.fill(HIST("hProtonCount"), 1.5);
            registry.fill(HIST("hDauProtonPt"), track.pt());
            registry.fill(HIST("hDauProtonMcPt"), mcparticle.pt());
            registry.fill(HIST("hDauProtonNsigmaProton"), track.tpcNSigmaPr());
            if (track.p() > 1.5) {
              registry.fill(HIST("hDauProtonNsigmaProtonAfterMomCut"), track.tpcNSigmaPr());
            }
            if (track.tpcNClsCrossedRows() < 70) {
              continue;
            }
            registry.fill(HIST("hProtonCount"), 2.5);
            if (TMath::Abs(track.eta()) > 0.9) {
              continue;
            }
            registry.fill(HIST("hProtonCount"), 3.5);
            if (track.pt() < minProtonPt || track.pt() > maxProtonPt) {
              continue;
            }
            registry.fill(HIST("hProtonCount"), 4.5);
            if (TMath::Abs(track.tpcNSigmaPr()) > 5) {
              continue;
            }
            registry.fill(HIST("hProtonCount"), 5.5);
          }
        }

        registry.fill(HIST("hProtonMcPt"), mcparticle.pt());
        registry.fill(HIST("hProtonMcP"), mcparticle.p());
        registry.fill(HIST("hProtonPt"), track.pt());
        registry.fill(HIST("hProtonP"), track.p());

        registry.fill(HIST("hProtonNsigmaProton"), track.tpcNSigmaPr());
        registry.fill(HIST("hProtonTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hProtonEta"), track.eta());
        registry.fill(HIST("hProtonMcRapidity"), mcparticle.y());
        registry.fill(HIST("hProtonTPCBB"), track.p() * track.sign(), track.tpcSignal());
      }

      // Pion
      if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
        registry.fill(HIST("hParticleCount"), 5.5);
        if (track.tpcNClsCrossedRows() > 70) {
          registry.fill(HIST("hParticleCount2"), 5.5);
          registry.fill(HIST("hPionTPCBBAfterTPCNclsCut"), track.p() * track.sign(), track.tpcSignal());
        }

        if (mcparticle.has_mothers()) {
          registry.fill(HIST("hPionCount"), 0.5);
          for (auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
            if (particleMother.pdgCode() != 1010010030 && particleMother.pdgCode() != -1010010030) {
              continue;
            }
            bool haveProton = false, havePion = false, haveDeuteron = false;
            bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
            for (auto& mcparticleDaughter : particleMother.daughters_as<aod::McParticles>()) {
              if (mcparticleDaughter.pdgCode() == 2212) {
                haveProton = true;
              }
              if (mcparticleDaughter.pdgCode() == -2212) {
                haveAntiProton = true;
              }
              if (mcparticleDaughter.pdgCode() == 211) {
                havePion = true;
              }
              if (mcparticleDaughter.pdgCode() == -211) {
                haveAntiPion = true;
              }
              if (mcparticleDaughter.pdgCode() == 1000010020) {
                haveDeuteron = true;
              }
              if (mcparticleDaughter.pdgCode() == -1000010020) {
                haveAntiDeuteron = true;
              }
            }
            if (!(haveProton && haveAntiPion && haveDeuteron && particleMother.pdgCode() == 1010010030) && !(haveAntiProton && havePion && haveAntiDeuteron && particleMother.pdgCode() == -1010010030)) {
              continue;
            }
            registry.fill(HIST("hPionCount"), 1.5);
            registry.fill(HIST("hDauPionPt"), track.pt());
            registry.fill(HIST("hDauPionMcPt"), mcparticle.pt());
            registry.fill(HIST("hDauPionNsigmaPion"), track.tpcNSigmaPi());
            if (track.p() > 0.2) {
              registry.fill(HIST("hDauPionNsigmaPionAfterMomCut"), track.tpcNSigmaPi());
            }
            if (track.tpcNClsCrossedRows() < 70) {
              continue;
            }
            registry.fill(HIST("hPionCount"), 2.5);
            if (TMath::Abs(track.eta()) > 0.9) {
              continue;
            }
            registry.fill(HIST("hPionCount"), 3.5);
            if (track.pt() < minPionPt || track.pt() > maxPionPt) {
              continue;
            }
            registry.fill(HIST("hPionCount"), 4.5);
            if (TMath::Abs(track.tpcNSigmaPi()) > 5) {
              continue;
            }
            registry.fill(HIST("hPionCount"), 5.5);
            if (TMath::Abs(track.dcaXY()) < dcapiontopv) {
              continue;
            }
            registry.fill(HIST("hPionCount"), 6.5);
          }
        }

        bool isFromHypertriton = false;
        if (mcparticle.has_mothers()) {
          for (auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
            if (particleMother.pdgCode() == 1010010030 || particleMother.pdgCode() == -1010010030) {
              isFromHypertriton = true;
            }
          }
        }

        while (!isFromHypertriton) {
          registry.fill(HIST("hBackgroundPionCount"), 0.5);
          if (track.tpcNClsCrossedRows() < 70) {
            break;
          }
          registry.fill(HIST("hBackgroundPionCount"), 1.5);
          if (TMath::Abs(track.eta()) > 0.9) {
            break;
          }
          registry.fill(HIST("hBackgroundPionCount"), 2.5);
          if (track.pt() < minPionPt || track.pt() > maxPionPt) {
            break;
          }
          registry.fill(HIST("hBackgroundPionCount"), 3.5);
          if (TMath::Abs(track.tpcNSigmaPi()) > 5) {
            break;
          }
          registry.fill(HIST("hBackgroundPionCount"), 4.5);
          break;
        }

        registry.fill(HIST("hPionMcPt"), mcparticle.pt());
        registry.fill(HIST("hPionMcP"), mcparticle.p());
        registry.fill(HIST("hPionPt"), track.pt());
        registry.fill(HIST("hPionP"), track.p());

        registry.fill(HIST("hPionNsigmaPion"), track.tpcNSigmaPi());
        registry.fill(HIST("hPionTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hPionEta"), track.eta());
        registry.fill(HIST("hPionMcRapidity"), mcparticle.y());
        registry.fill(HIST("hPionTPCBB"), track.p() * track.sign(), track.tpcSignal());
      }

      // Deuteron
      if (mcparticle.pdgCode() == 1000010020 || mcparticle.pdgCode() == -1000010020) {
        registry.fill(HIST("hParticleCount"), 6.5);
        if (track.tpcNClsCrossedRows() > 70) {
          registry.fill(HIST("hParticleCount2"), 6.5);
          registry.fill(HIST("hDeuteronTPCBBAfterTPCNclsCut"), track.p() * track.sign(), track.tpcSignal());
        }

        if (mcparticle.has_mothers()) {
          registry.fill(HIST("hDeuteronCount"), 0.5);
          for (auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
            if (particleMother.pdgCode() != 1010010030 && particleMother.pdgCode() != -1010010030) {
              continue;
            }
            bool haveProton = false, havePion = false, haveDeuteron = false;
            bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
            for (auto& mcparticleDaughter : particleMother.daughters_as<aod::McParticles>()) {
              if (mcparticleDaughter.pdgCode() == 2212) {
                haveProton = true;
              }
              if (mcparticleDaughter.pdgCode() == -2212) {
                haveAntiProton = true;
              }
              if (mcparticleDaughter.pdgCode() == 211) {
                havePion = true;
              }
              if (mcparticleDaughter.pdgCode() == -211) {
                haveAntiPion = true;
              }
              if (mcparticleDaughter.pdgCode() == 1000010020) {
                haveDeuteron = true;
              }
              if (mcparticleDaughter.pdgCode() == -1000010020) {
                haveAntiDeuteron = true;
              }
            }
            if (!(haveProton && haveAntiPion && haveDeuteron && particleMother.pdgCode() == 1010010030) && !(haveAntiProton && havePion && haveAntiDeuteron && particleMother.pdgCode() == -1010010030)) {
              continue;
            }
            registry.fill(HIST("hDeuteronCount"), 1.5);
            registry.fill(HIST("hDauDeuteronPt"), track.pt());
            registry.fill(HIST("hDauDeuteronMcPt"), mcparticle.pt());
            registry.fill(HIST("hDauDeuteronNsigmaDeuteron"), track.tpcNSigmaDe());
            if (track.p() > 1.5) {
              registry.fill(HIST("hDauDeuteronNsigmaDeuteronAfterMomCut"), track.tpcNSigmaDe());
            }
            if (track.tpcNClsCrossedRows() < 70) {
              continue;
            }
            registry.fill(HIST("hDeuteronCount"), 2.5);
            if (TMath::Abs(track.eta()) > 0.9) {
              continue;
            }
            registry.fill(HIST("hDeuteronCount"), 3.5);
            if (track.pt() < minDeuteronPt || track.pt() > maxDeuteronPt) {
              continue;
            }
            registry.fill(HIST("hDeuteronCount"), 4.5);
            if (TMath::Abs(track.tpcNSigmaDe()) > 5) {
              continue;
            }
            registry.fill(HIST("hDeuteronCount"), 5.5);
          }
        }

        registry.fill(HIST("hDeuteronMcPt"), mcparticle.pt());
        registry.fill(HIST("hDeuteronMcP"), mcparticle.p());
        registry.fill(HIST("hDeuteronPt"), track.pt());
        registry.fill(HIST("hDeuteronP"), track.p());

        registry.fill(HIST("hDeuteronNsigmaDeuteron"), track.tpcNSigmaDe());
        registry.fill(HIST("hDeuteronTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hDeuteronEta"), track.eta());
        registry.fill(HIST("hDeuteronMcRapidity"), mcparticle.y());
        registry.fill(HIST("hDeuteronTPCBB"), track.p() * track.sign(), track.tpcSignal());
      }
    }
  }
};

namespace o2::aod
{
namespace mcparticleidcheck
{
DECLARE_SOA_COLUMN(HypertritonDaughterPdgCode, hypertritonDaughterPdgCode, int);
}
DECLARE_SOA_TABLE(McHypertritonCheck, "AOD", "HypCheck",
                  mcparticleidcheck::HypertritonDaughterPdgCode)
} // namespace o2::aod

// check the performance of mcparticle
struct hypertriton3bodyMcParticleCount {
  // Basic checks
  Produces<aod::McHypertritonCheck> hypertrtitonCheckTable;
  HistogramRegistry registry{
    "registry",
    {
      {"hTotalMcCollCounter", "hTotalMcCollCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},

      {"h3dMCHypertriton", "h3dMCHypertriton", {HistType::kTH3F, {{20, -1.0f, 1.0f, "Rapidity"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {50, 0.0f, 50.0f, "ct(cm)"}}}},
      {"h3dMCDecayedHypertriton", "h3dMCDecayedHypertriton", {HistType::kTH3F, {{20, -1.0f, 1.0f, "Rapidity"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {50, 0.0f, 50.0f, "ct(cm)"}}}},
      {"hMcHypertritonCheck", "hMcHypertritonCheck", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
      {"hMcPhysicalPrimaryParticleCount", "hMcPhysicalPrimaryParticleCount", {HistType::kTH1F, {{8, 0.0f, 8.0f}}}},

      {"hMcHypertritonCount", "hMcHypertritonCount", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hMcHypertritonDecayCount", "hMcHypertritonDecayCount", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hMcHypertritonPt", "hMcHypertritonPt", {HistType::kTH1F, {{300, 0.0f, 15.0f}}}},
      {"hMcHypertritonLifeTime", "hMcHypertritonLifeTime", {HistType::kTH1F, {{500, 0.0f, 500.0f}}}},
      {"hMcProtonPt", "hMcProtonPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hMcPionPt", "hMcPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hMcDeuteronPt", "hMcDeuteronPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},

      {"hMcRecoInvMass", "hMcRecoInvMass", {HistType::kTH1F, {{100, 2.95, 3.05f}}}},
    },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hTotalMcCollCounter"))->GetXaxis()->SetBinLabel(1, "Total Count");
    registry.get<TH1>(HIST("hTotalMcCollCounter"))->GetXaxis()->SetBinLabel(2, "Recoonstructed");

    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(1, "PdgCode");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(2, "Rapidity");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(3, "Lifetime");
    registry.get<TH1>(HIST("hMcHypertritonCheck"))->GetXaxis()->SetBinLabel(4, "PtCut");

    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(1, "ReadinAfterEvtRecSelCut(off)");
    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(2, "IsPhysicalPrimary");
    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(3, "y<0.9(off)");
    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(4, "(Anti)Proton");
    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(5, "(Anti)Pion");
    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(6, "(Anti)Deuteron");
    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(7, "(Anti)Hypertriton");
    registry.get<TH1>(HIST("hMcPhysicalPrimaryParticleCount"))->GetXaxis()->SetBinLabel(8, "HasDaughter");
    registry.get<TH1>(HIST("hMcHypertritonCount"))->GetXaxis()->SetBinLabel(1, "Generated number");
    registry.get<TH1>(HIST("hMcHypertritonCount"))->GetXaxis()->SetBinLabel(2, "primary Hypertriton mothers");
    registry.get<TH1>(HIST("hMcHypertritonDecayCount"))->GetXaxis()->SetBinLabel(1, "confirm to 3-body decay");
    registry.get<TH1>(HIST("hMcHypertritonDecayCount"))->GetXaxis()->SetBinLabel(2, "Hypertriton Candidates");
    registry.get<TH1>(HIST("hMcHypertritonDecayCount"))->GetXaxis()->SetBinLabel(3, "AntiHypertriton Candidates");
  }

  Configurable<float> rapidityMCcut{"rapidityMCcut", 0.9, "rapidity cut MC count"};
  Configurable<bool> eventSelectionMC{"eventSelectionMC", true, "event selection MC count"};

  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>& collisions)
  {
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (eventSelectionMC && !collision.sel8()) {
        continue;
      }
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);

    registry.fill(HIST("hTotalMcCollCounter"), 0.5);

    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();
    if (evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      registry.fill(HIST("hTotalMcCollCounter"), 1.5);
      return;
    }

    for (auto& mcparticle : mcParticles) {

      if (mcparticle.pdgCode() == 2212 || mcparticle.pdgCode() == -2212) {
        registry.fill(HIST("hMcProtonPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
        registry.fill(HIST("hMcPionPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == 1000010020 || mcparticle.pdgCode() == -1000010020) {
        registry.fill(HIST("hMcDeuteronPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hMcHypertritonCount"), 0.5);
        registry.fill(HIST("hMcHypertritonCheck"), 0.5);
        registry.fill(HIST("hMcHypertritonPt"), mcparticle.pt());

        // double dauProtonPos[3] = {-999, -999, -999};
        // double dauPionPos[3] = {-999, -999, -999};
        double dauDeuteronPos[3] = {-999, -999, -999};
        double dauProtonMom[3] = {-999, -999, -999};
        double dauPionMom[3] = {-999, -999, -999};
        double dauDeuteronMom[3] = {-999, -999, -999};
        bool haveProton = false, havePion = false, haveDeuteron = false;
        bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
        for (auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          hypertrtitonCheckTable(mcparticleDaughter.pdgCode());
          if (mcparticleDaughter.pdgCode() == 2212) {
            haveProton = true;
          }
          if (mcparticleDaughter.pdgCode() == -2212) {
            haveAntiProton = true;
          }
          if (mcparticleDaughter.pdgCode() == 211) {
            havePion = true;
          }
          if (mcparticleDaughter.pdgCode() == -211) {
            haveAntiPion = true;
          }
          if (mcparticleDaughter.pdgCode() == 1000010020) {
            haveDeuteron = true;
          }
          if (mcparticleDaughter.pdgCode() == -1000010020) {
            haveAntiDeuteron = true;
          }
          if (std::abs(mcparticleDaughter.pdgCode()) == 2212) {
            /*dauProtonPos[0] = mcparticleDaughter.vx();
              dauProtonPos[1] = mcparticleDaughter.vy();
              dauProtonPos[2] = mcparticleDaughter.vz();*/
            dauProtonMom[0] = mcparticleDaughter.px();
            dauProtonMom[1] = mcparticleDaughter.py();
            dauProtonMom[2] = mcparticleDaughter.pz();
          }
          if (std::abs(mcparticleDaughter.pdgCode()) == 211) {
            /*dauPionPos[0] = mcparticleDaughter.vx();
              dauPionPos[1] = mcparticleDaughter.vy();
              dauPionPos[2] = mcparticleDaughter.vz();*/
            dauPionMom[0] = mcparticleDaughter.px();
            dauPionMom[1] = mcparticleDaughter.py();
            dauPionMom[2] = mcparticleDaughter.pz();
          }
          if (std::abs(mcparticleDaughter.pdgCode()) == 1000010020) {
            dauDeuteronPos[0] = mcparticleDaughter.vx();
            dauDeuteronPos[1] = mcparticleDaughter.vy();
            dauDeuteronPos[2] = mcparticleDaughter.vz();
            dauDeuteronMom[0] = mcparticleDaughter.px();
            dauDeuteronMom[1] = mcparticleDaughter.py();
            dauDeuteronMom[2] = mcparticleDaughter.pz();
          }
        }
        if (haveProton && haveAntiPion && haveDeuteron && mcparticle.pdgCode() == 1010010030) {
          registry.fill(HIST("hMcHypertritonDecayCount"), 0.5);
          registry.fill(HIST("hMcHypertritonDecayCount"), 1.5);
        }
        if (haveAntiProton && havePion && haveAntiDeuteron && mcparticle.pdgCode() == -1010010030) {
          registry.fill(HIST("hMcHypertritonDecayCount"), 0.5);
          registry.fill(HIST("hMcHypertritonDecayCount"), 2.5);
        }
        double MClifetime = RecoDecay::sqrtSumOfSquares(dauDeuteronPos[0] - mcparticle.vx(), dauDeuteronPos[1] - mcparticle.vy(), dauDeuteronPos[2] - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
        if ((haveProton && haveAntiPion && haveDeuteron && mcparticle.pdgCode() == 1010010030) || (haveAntiProton && havePion && haveAntiDeuteron && mcparticle.pdgCode() == -1010010030)) {
          registry.fill(HIST("hMcRecoInvMass"), RecoDecay::m(array{array{dauProtonMom[0], dauProtonMom[1], dauProtonMom[2]}, array{dauPionMom[0], dauPionMom[1], dauPionMom[2]}, array{dauDeuteronMom[0], dauDeuteronMom[1], dauDeuteronMom[2]}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron}));
          registry.fill(HIST("h3dMCDecayedHypertriton"), mcparticle.y(), mcparticle.pt(), MClifetime);
        }

        int daughterPionCount = 0;
        for (auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          if (std::abs(mcparticleDaughter.pdgCode()) == 211) {
            daughterPionCount++;
          }
        }

        registry.fill(HIST("hMcHypertritonLifeTime"), MClifetime);
        registry.fill(HIST("h3dMCHypertriton"), mcparticle.y(), mcparticle.pt(), MClifetime);

        // Count for hypertriton N_gen
        if (TMath::Abs(mcparticle.y()) < 1) {
          registry.fill(HIST("hMcHypertritonCheck"), 1.5);
          if (MClifetime < 40) {
            registry.fill(HIST("hMcHypertritonCheck"), 2.5);
            if (mcparticle.pt() > 1 && mcparticle.pt() < 10) {
              registry.fill(HIST("hMcHypertritonCheck"), 3.5);
            }
          }
        }
      }

      registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 0.5);
      if (!mcparticle.isPhysicalPrimary()) {
        continue;
      }
      registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 1.5);
      // if (TMath::Abs(mcparticle.y()) > rapidityMCcut) {continue;}
      registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 2.5);

      if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
        registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 3.5);
      } else if (mcparticle.pdgCode() == 2212 || mcparticle.pdgCode() == -2212) {
        registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 4.5);
      } else if (mcparticle.pdgCode() == 1000010020 || mcparticle.pdgCode() == -1000010020) {
        registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 5.5);
      } else if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
        registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 6.5);
        registry.fill(HIST("hMcHypertritonCount"), 1.5);
      }

      if (!mcparticle.has_daughters()) {
        continue;
      }
      registry.fill(HIST("hMcPhysicalPrimaryParticleCount"), 7.5);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertriton3bodyAnalysisMc>(cfgc),
    adaptAnalysisTask<hypertriton3bodyLabelCheck>(cfgc),
    adaptAnalysisTask<hypertriton3bodyTrackMcinfo>(cfgc),
    adaptAnalysisTask<hypertriton3bodyMcParticleCount>(cfgc),
    adaptAnalysisTask<hypertriton3bodyQaMc>(cfgc),
  };
}
