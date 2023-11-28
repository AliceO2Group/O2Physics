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
//

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

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::pidTPCFullPr, aod::pidTOFFullDe>;

struct hypertriton3bodyQa {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxRadius", "hVtxRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hVtxCosPA", "hVtxCosPA", {HistType::kTH1F, {{1000, 0.9f, 1.0f}}}},
      {"hPtProton", "hPtProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiPion", "hPtAntiPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtDeuteron", "hPtDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiProton", "hPtAntiProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtPion", "hPtPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiDeuteron", "hPtAntiDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDCAProtonToPV", "hDCAProtonToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAPionToPV", "hDCAPionToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCADeuteronToPV", "hDCADeuteronToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hProtonTPCNcls", "hProtonTPCNcls", {HistType::kTH1F, {{300, 0, 300, "TPC cluster"}}}},
      {"hPionTPCNcls", "hPionTPCNcls", {HistType::kTH1F, {{300, 0, 300, "TPC cluster"}}}},
      {"hDeuteronTPCNcls", "hDeuteronTPCNcls", {HistType::kTH1F, {{300, 0, 300, "TPC cluster"}}}},
      {"hDCAVtxDau", "hDCAVtxDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hVtxPt", "hVtxPt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTrack0Pt", "hTrack0Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTrack1Pt", "hTrack1Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTrack2Pt", "hTrack2Pt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTOFPIDDeuteron", "hTOFPIDDeuteron", {HistType::kTH1F, {{2000, -100.0f, 100.0f}}}},
      {"hDeuTOFNsigma", "Deuteron TOF Nsigma distribution", {HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {2000, -100, 100, "TOF n#sigma"}}}},
      {"hTOFPIDDeuteronWithTPC", "hTOFPIDDeuteronWithTPC", {HistType::kTH1F, {{2000, -100.0f, 100.0f}}}},
      {"hDeuTOFNsigmaWithTPC", "Deuteron TOF Nsigma distribution", {HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {2000, -100, 100, "TOF n#sigma"}}}},
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
      auto track0 = vtx.track0_as<MyTracks>();
      auto track1 = vtx.track1_as<MyTracks>();
      auto track2 = vtx.track2_as<MyTracks>();

      registry.fill(HIST("hVtxRadius"), vtx.vtxradius());
      registry.fill(HIST("hVtxCosPA"), vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAVtxDau"), vtx.dcaVtxdaughters());
      registry.fill(HIST("hVtxPt"), vtx.pt());
      registry.fill(HIST("hTrack0Pt"), vtx.track0pt());
      registry.fill(HIST("hTrack1Pt"), vtx.track1pt());
      registry.fill(HIST("hTrack2Pt"), vtx.track2pt());
      registry.fill(HIST("hMassHypertriton"), vtx.mHypertriton());
      registry.fill(HIST("hMassAntiHypertriton"), vtx.mAntiHypertriton());
      registry.fill(HIST("hTOFPIDDeuteron"), track2.tofNSigmaDe());
      registry.fill(HIST("hDeuTOFNsigma"), track2.tpcInnerParam() * track2.sign(), track2.tofNSigmaDe());
      if (std::abs(track2.tpcNSigmaDe()) < 5) {
        registry.fill(HIST("hTOFPIDDeuteronWithTPC"), track2.tofNSigmaDe());
        registry.fill(HIST("hDeuTOFNsigmaWithTPC"), track2.tpcInnerParam() * track2.sign(), track2.tofNSigmaDe());
      }
      if (track2.sign() > 0) {
        registry.fill(HIST("hPtProton"), track0.pt());
        registry.fill(HIST("hPtAntiPion"), track1.pt());
        registry.fill(HIST("hPtDeuteron"), track2.pt());
        registry.fill(HIST("hDCAProtonToPV"), vtx.dcatrack0topv());
        registry.fill(HIST("hDCAPionToPV"), vtx.dcatrack1topv());
        registry.fill(HIST("hProtonTPCNcls"), track0.tpcNClsCrossedRows());
        registry.fill(HIST("hPionTPCNcls"), track1.tpcNClsCrossedRows());
      } else {
        registry.fill(HIST("hPtPion"), track0.pt());
        registry.fill(HIST("hPtAntiProton"), track1.pt());
        registry.fill(HIST("hPtAntiDeuteron"), track2.pt());
        registry.fill(HIST("hDCAProtonToPV"), vtx.dcatrack1topv());
        registry.fill(HIST("hDCAPionToPV"), vtx.dcatrack0topv());
        registry.fill(HIST("hProtonTPCNcls"), track1.tpcNClsCrossedRows());
        registry.fill(HIST("hPionTPCNcls"), track0.tpcNClsCrossedRows());
      }
      registry.fill(HIST("hDCADeuteronToPV"), vtx.dcatrack2topv());
      registry.fill(HIST("hDeuteronTPCNcls"), track2.tpcNClsCrossedRows());
    }
  }
};

struct hypertriton3bodyAnalysis {

  // Selection criteria
  Configurable<double> vtxcospa{"vtxcospa", 0.9, "Vtx CosPA"};          // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"}; // loose cut
  Configurable<float> dcapiontopv{"dcapiontopv", .00, "DCA Pion To PV"};
  Configurable<float> etacut{"etacut", 1, "etacut"};
  Configurable<float> rapidity{"rapidity", 1, "rapidity"};
  Configurable<float> TofPidNsigmaMin{"TofPidNsigmaMin", -4, "TofPidNsigmaMin"};
  Configurable<float> TofPidNsigmaMax{"TofPidNsigmaMax", 8, "TofPidNsigmaMax"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<float> lifetimecut{"lifetimecut", 40., "lifetimecut"}; // ct
  Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
  Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
  Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
  Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
  Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
  Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};
  Configurable<float> minDeuteronPUseTOF{"minDeuteronPUseTOF", 1, "minDeuteronPt Enable TOF PID"};
  Configurable<float> h3LMassLowerlimit{"h3LMassLowerlimit", 2.96, "Hypertriton mass lower limit"};
  Configurable<float> h3LMassUpperlimit{"h3LMassUpperlimit", 3.04, "Hypertriton mass upper limit"};
  Configurable<int> mincrossedrowsproton{"mincrossedrowsproton", 90, "min tpc crossed rows for pion"};
  Configurable<int> mincrossedrowspion{"mincrossedrowspion", 70, "min tpc crossed rows"};
  Configurable<int> mincrossedrowsdeuteron{"mincrossedrowsdeuteron", 100, "min tpc crossed rows for deuteron"};

  Configurable<float> mcsigma{"mcsigma", 0.0015, "sigma of mc invariant mass fit"}; // obtained from MC

  HistogramRegistry registry{
    "registry",
    {
      {"hSelectedEventCounter", "hSelectedEventCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hSelectedCandidatesCounter", "hSelectedCandidatesCounter", {HistType::kTH1F, {{11, 0.0f, 11.0f}}}},
      {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
      {"hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
      {"hMassHypertritonTotal", "hMassHypertritonTotal", {HistType::kTH1F, {{300, 2.9f, 3.2f}}}},
      {"hPtProton", "hPtProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiPion", "hPtAntiPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtDeuteron", "hPtDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiProton", "hPtAntiProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtPion", "hPtPion", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiDeuteron", "hPtAntiDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDCAProtonToPV", "hDCAProtonToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAPionToPV", "hDCAPionToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCADeuteronToPV", "hDCADeuteronToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hProtonTPCNcls", "hProtonTPCNcls", {HistType::kTH1F, {{180, 0, 180, "TPC cluster"}}}},
      {"hPionTPCNcls", "hPionTPCNcls", {HistType::kTH1F, {{180, 0, 180, "TPC cluster"}}}},
      {"hDeuteronTPCNcls", "hDeuteronTPCNcls", {HistType::kTH1F, {{180, 0, 180, "TPC cluster"}}}},
      {"hVtxCosPA", "hVtxCosPA", {HistType::kTH1F, {{1000, 0.9f, 1.0f}}}},
      {"hDCAVtxDau", "hDCAVtxDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hTOFPIDDeuteron", "hTOFPIDDeuteron", {HistType::kTH1F, {{2000, -100.0f, 100.0f}}}},
      {"hTPCPIDProton", "hTPCPIDProton", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hTPCPIDPion", "hTPCPIDPion", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hTPCPIDDeuteron", "hTPCPIDDeuteron", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
      {"hProtonTPCBB", "hProtonTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hPionTPCBB", "hPionTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDeuteronTPCBB", "hDeuteronTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hProtonTPCVsPt", "hProtonTPCVsPt", {HistType::kTH2F, {{50, 0.0f, 5.0f, "#it{p}_{T} (GeV/c)"}, {240, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hPionTPCVsPt", "hPionTPCVsPt", {HistType::kTH2F, {{20, 0.0f, 2.0f, "#it{p}_{T} (GeV/c)"}, {240, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDeuteronTPCVsPt", "hDeuteronTPCVsPt", {HistType::kTH2F, {{80, 0.0f, 8.0f, "#it{p}_{T} (GeV/c)"}, {240, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDalitz", "hDalitz", {HistType::kTH2F, {{120, 7.85, 8.45, "M^{2}(dp) (GeV^{2}/c^{4})"}, {60, 1.1, 1.4, "M^{2}(p#pi) (GeV^{2}/c^{4})"}}}},
      {"h3dMassHypertriton", "h3dMassHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dMassAntiHypertriton", "h3dMassAntiHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dTotalHypertriton", "h3dTotalHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
    },
  };

  //------------------------------------------------------------------
  // Fill stats histograms
  enum vtxstep { kCandAll = 0,
                 kCandCosPA,
                 kCandDauEta,
                 kCandRapidity,
                 kCandct,
                 kCandDcaDau,
                 kCandTOFPID,
                 kCandTPCPID,
                 kCandTPCNcls,
                 kCandDauPt,
                 kCandDcaToPV,
                 kNVtxSteps };

  Configurable<int> saveDcaHist{"saveDcaHist", 1, "saveDcaHist"};
  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};

  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisHypertriton = {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"};

    if (saveDcaHist == 1) {
      registry.add("h3dMassHypertritonDca", "h3dMassHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
      registry.add("h3dMassAntiHypertritonDca", "h3dMassAntiHypertritonDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisHypertriton}});
    }

    TString CandCounterbinLabel[11] = {"Total", "VtxCosPA", "TrackEta", "MomRapidity", "Lifetime", "DcaV0Dau", "d TOFPID", "TPCPID&Mass", "TPCNcls", "DauPt", "PionDcatoPV"};
    for (int i{0}; i < kNVtxSteps; i++) {
      registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(i + 1, CandCounterbinLabel[i]);
    }
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, MyTracks const& tracks)
  {
    registry.fill(HIST("hSelectedEventCounter"), 0.5);
    if (eventSelection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hSelectedEventCounter"), 1.5);

    bool if_hasvtx = false;

    for (auto& vtx : vtx3bodydatas) {

      auto track0 = vtx.track0_as<MyTracks>();
      auto track1 = vtx.track1_as<MyTracks>();
      auto track2 = vtx.track2_as<MyTracks>();

      registry.fill(HIST("hSelectedCandidatesCounter"), kCandAll);
      if (vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()) < vtxcospa) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), kCandCosPA);
      if (TMath::Abs(track0.eta()) > etacut || TMath::Abs(track1.eta()) > etacut || TMath::Abs(track2.eta()) > etacut) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), kCandDauEta);
      if (TMath::Abs(vtx.yHypertriton()) > rapidity) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), kCandRapidity);
      double ct = vtx.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassHyperTriton;
      if (ct > lifetimecut) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), kCandct);
      if (vtx.dcaVtxdaughters() > dcavtxdau) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), kCandDcaDau);
      if ((track2.tofNSigmaDe() < TofPidNsigmaMin || track2.tofNSigmaDe() > TofPidNsigmaMax) && track2.p() > minDeuteronPUseTOF) {
        continue;
      }
      registry.fill(HIST("hSelectedCandidatesCounter"), kCandTOFPID);

      // 3sigma region for Dalitz plot
      double lowersignallimit = o2::constants::physics::MassHyperTriton - 3 * mcsigma;
      double uppersignallimit = o2::constants::physics::MassHyperTriton + 3 * mcsigma;

      // Hypertriton
      if (TMath::Abs(track0.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(track1.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(track2.tpcNSigmaDe()) < TpcPidNsigmaCut && vtx.mHypertriton() > h3LMassLowerlimit && vtx.mHypertriton() < h3LMassUpperlimit) {

        registry.fill(HIST("hSelectedCandidatesCounter"), kCandTPCPID);

        if (track0.tpcNClsCrossedRows() > mincrossedrowsproton && track1.tpcNClsCrossedRows() > mincrossedrowspion && track2.tpcNClsCrossedRows() > mincrossedrowsdeuteron) {

          registry.fill(HIST("hSelectedCandidatesCounter"), kCandTPCNcls);

          if (vtx.track0pt() > minProtonPt && vtx.track0pt() < maxProtonPt && vtx.track1pt() > minPionPt && vtx.track1pt() < maxPionPt && vtx.track2pt() > minDeuteronPt && vtx.track2pt() < maxDeuteronPt) {
            registry.fill(HIST("hSelectedCandidatesCounter"), kCandDauPt);

            if (TMath::Abs(vtx.dcatrack1topv()) > dcapiontopv) {
              if_hasvtx = true;
              registry.fill(HIST("hSelectedCandidatesCounter"), kCandDcaToPV);

              registry.fill(HIST("hVtxCosPA"), vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()));
              registry.fill(HIST("hDCAVtxDau"), vtx.dcaVtxdaughters());
              registry.fill(HIST("hPtProton"), vtx.track0pt());
              registry.fill(HIST("hPtAntiPion"), vtx.track1pt());
              registry.fill(HIST("hPtDeuteron"), vtx.track2pt());
              registry.fill(HIST("hDCAProtonToPV"), vtx.dcatrack0topv());
              registry.fill(HIST("hDCAPionToPV"), vtx.dcatrack1topv());
              registry.fill(HIST("hDCADeuteronToPV"), vtx.dcatrack2topv());
              registry.fill(HIST("hProtonTPCNcls"), track0.tpcNClsCrossedRows());
              registry.fill(HIST("hPionTPCNcls"), track1.tpcNClsCrossedRows());
              registry.fill(HIST("hDeuteronTPCNcls"), track2.tpcNClsCrossedRows());
              registry.fill(HIST("hTOFPIDDeuteron"), track2.tofNSigmaDe());
              registry.fill(HIST("hTPCPIDProton"), track0.tpcNSigmaPr());
              registry.fill(HIST("hTPCPIDPion"), track1.tpcNSigmaPi());
              registry.fill(HIST("hTPCPIDDeuteron"), track2.tpcNSigmaDe());
              registry.fill(HIST("hProtonTPCBB"), track0.p(), track0.tpcSignal());
              registry.fill(HIST("hPionTPCBB"), -track1.p(), track1.tpcSignal());
              registry.fill(HIST("hDeuteronTPCBB"), track2.p(), track0.tpcSignal());
              registry.fill(HIST("hProtonTPCVsPt"), vtx.track0pt(), track0.tpcNSigmaPr());
              registry.fill(HIST("hPionTPCVsPt"), vtx.track1pt(), track1.tpcNSigmaPi());
              registry.fill(HIST("hDeuteronTPCVsPt"), vtx.track2pt(), track2.tpcNSigmaDe());
              registry.fill(HIST("hMassHypertriton"), vtx.mHypertriton());
              registry.fill(HIST("hMassHypertritonTotal"), vtx.mHypertriton());
              registry.fill(HIST("h3dMassHypertriton"), 0., vtx.pt(), vtx.mHypertriton()); // collision.centV0M() instead of 0. once available
              registry.fill(HIST("h3dTotalHypertriton"), ct, vtx.pt(), vtx.mHypertriton());
              if (vtx.mHypertriton() > lowersignallimit && vtx.mHypertriton() < uppersignallimit) {
                registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{vtx.pxtrack0(), vtx.pytrack0(), vtx.pztrack0()}, array{vtx.pxtrack2(), vtx.pytrack2(), vtx.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}),
                              RecoDecay::m2(array{array{vtx.pxtrack0(), vtx.pytrack0(), vtx.pztrack0()}, array{vtx.pxtrack1(), vtx.pytrack1(), vtx.pztrack1()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
              }

              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassHypertritonDca"), vtx.dcaVtxdaughters(), vtx.pt(), vtx.mHypertriton());
              }
            }
          }
        }
      }

      // AntiHypertriton
      if (TMath::Abs(track0.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(track1.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(track2.tpcNSigmaDe()) < TpcPidNsigmaCut && vtx.mAntiHypertriton() > h3LMassLowerlimit && vtx.mAntiHypertriton() < h3LMassUpperlimit) {

        registry.fill(HIST("hSelectedCandidatesCounter"), kCandTPCPID);

        if (track0.tpcNClsCrossedRows() > mincrossedrowspion && track1.tpcNClsCrossedRows() > mincrossedrowsproton && track2.tpcNClsCrossedRows() > mincrossedrowsdeuteron) {

          registry.fill(HIST("hSelectedCandidatesCounter"), kCandTPCNcls);

          if (vtx.track0pt() > minPionPt && vtx.track0pt() < maxPionPt && vtx.track1pt() > minProtonPt && vtx.track1pt() < maxProtonPt && vtx.track2pt() > minDeuteronPt && vtx.track2pt() < maxDeuteronPt) {
            registry.fill(HIST("hSelectedCandidatesCounter"), kCandDauPt);
            if (TMath::Abs(vtx.dcatrack0topv()) > dcapiontopv) {
              if_hasvtx = true;
              registry.fill(HIST("hSelectedCandidatesCounter"), kCandDcaToPV);

              registry.fill(HIST("hVtxCosPA"), vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()));
              registry.fill(HIST("hDCAVtxDau"), vtx.dcaVtxdaughters());
              registry.fill(HIST("hPtAntiProton"), vtx.track1pt());
              registry.fill(HIST("hPtPion"), vtx.track0pt());
              registry.fill(HIST("hPtAntiDeuteron"), vtx.track2pt());
              registry.fill(HIST("hDCAProtonToPV"), vtx.dcatrack1topv());
              registry.fill(HIST("hDCAPionToPV"), vtx.dcatrack0topv());
              registry.fill(HIST("hDCADeuteronToPV"), vtx.dcatrack2topv());
              registry.fill(HIST("hProtonTPCNcls"), track1.tpcNClsCrossedRows());
              registry.fill(HIST("hPionTPCNcls"), track0.tpcNClsCrossedRows());
              registry.fill(HIST("hDeuteronTPCNcls"), track2.tpcNClsCrossedRows());
              registry.fill(HIST("hTOFPIDDeuteron"), track2.tofNSigmaDe());
              registry.fill(HIST("hTPCPIDProton"), track1.tpcNSigmaPr());
              registry.fill(HIST("hTPCPIDPion"), track0.tpcNSigmaPi());
              registry.fill(HIST("hTPCPIDDeuteron"), track2.tpcNSigmaDe());
              registry.fill(HIST("hProtonTPCBB"), -track1.p(), track1.tpcSignal());
              registry.fill(HIST("hPionTPCBB"), track0.p(), track0.tpcSignal());
              registry.fill(HIST("hDeuteronTPCBB"), -track2.p(), track0.tpcSignal());
              registry.fill(HIST("hProtonTPCVsPt"), vtx.track1pt(), track1.tpcNSigmaPr());
              registry.fill(HIST("hPionTPCVsPt"), vtx.track0pt(), track0.tpcNSigmaPi());
              registry.fill(HIST("hDeuteronTPCVsPt"), vtx.track2pt(), track2.tpcNSigmaDe());
              registry.fill(HIST("hMassAntiHypertriton"), vtx.mAntiHypertriton());
              registry.fill(HIST("hMassHypertritonTotal"), vtx.mAntiHypertriton());
              registry.fill(HIST("h3dMassAntiHypertriton"), 0., vtx.pt(), vtx.mAntiHypertriton()); // collision.centV0M() instead of 0. once available
              registry.fill(HIST("h3dTotalHypertriton"), ct, vtx.pt(), vtx.mAntiHypertriton());
              if (vtx.mAntiHypertriton() > lowersignallimit && vtx.mAntiHypertriton() < uppersignallimit) {
                registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{vtx.pxtrack1(), vtx.pytrack1(), vtx.pztrack1()}, array{vtx.pxtrack2(), vtx.pytrack2(), vtx.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}),
                              RecoDecay::m2(array{array{vtx.pxtrack1(), vtx.pytrack1(), vtx.pztrack1()}, array{vtx.pxtrack0(), vtx.pytrack0(), vtx.pztrack0()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
              }

              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassAntiHypertritonDca"), vtx.dcaVtxdaughters(), vtx.pt(), vtx.mAntiHypertriton());
              }
            }
          }
        }
      }
    }

    if (if_hasvtx)
      registry.fill(HIST("hSelectedEventCounter"), 2.5);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertriton3bodyAnalysis>(cfgc),
    adaptAnalysisTask<hypertriton3bodyQa>(cfgc),
  };
}
