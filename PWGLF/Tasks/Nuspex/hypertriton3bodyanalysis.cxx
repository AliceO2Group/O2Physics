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
// StoredVtx3BodyDatas analysis task
// ========================
//
// This code loops over a StoredVtx3BodyDatas table and produces some
// standard analysis output. It requires either
// the hypertriton3bodybuilder or hypertriton3bodyfinder (not recommended) tasks
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

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::pidTOFFullDe>;
using MCLabeledFullTracksExtIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

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
      {"hDeuTOFNsigmaWithTPC", "Deuteron TOF Nsigma distribution", {HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {1000, -100, 100, "TOF n#sigma"}}}},
    },
  };
  void init(InitContext const&)
  {
    AxisSpec massAxis = {120, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"};
    registry.add("hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {massAxis}});
    registry.add("hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {massAxis}});
  }
  void process(aod::Collision const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, FullTracksExtIU const& tracks)
  {
    for (auto& vtx : vtx3bodydatas) {
      auto track0 = vtx.track0_as<FullTracksExtIU>();
      auto track1 = vtx.track1_as<FullTracksExtIU>();
      auto track2 = vtx.track2_as<FullTracksExtIU>();

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
  Configurable<double> vtxcospa{"vtxcospa", 0.99, "Vtx CosPA"};          // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"}; // loose cut
  Configurable<float> dcapiontopv{"dcapiontopv", .00, "DCA Pion To PV"};
  Configurable<float> etacut{"etacut", 1, "etacut"};
  Configurable<float> rapidity{"rapidity", 1, "rapidity"};
  Configurable<float> TofPidNsigmaMin{"TofPidNsigmaMin", -4, "TofPidNsigmaMin"};
  Configurable<float> TofPidNsigmaMax{"TofPidNsigmaMax", 4, "TofPidNsigmaMax"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};
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
      {"hDeuteronTOFVsPBeforeTOFCut", "hDeuteronTOFVsPBeforeTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronTOFVsPAtferTOFCut", "hDeuteronTOFVsPAtferTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},

      {"hDalitz", "hDalitz", {HistType::kTH2F, {{120, 7.85, 8.45, "M^{2}(dp) (GeV^{2}/c^{4})"}, {60, 1.1, 1.4, "M^{2}(p#pi) (GeV^{2}/c^{4})"}}}},
      {"h3dMassHypertriton", "h3dMassHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dMassAntiHypertriton", "h3dMassAntiHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dTotalHypertriton", "h3dTotalHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},

      {"hSelectedTrueHypertritonCounter", "hSelectedTrueHypertritonCounter", {HistType::kTH1F, {{11, 0.0f, 11.0f}}}},
      {"hDeuteronTOFVsPBeforeTOFCutSig", "hDeuteronTOFVsPBeforeTOFCutSig", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronTOFVsPAtferTOFCutSig", "hDeuteronTOFVsPAtferTOFCutSig", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"h3dTotalTrueHypertriton", "h3dTotalTrueHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
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
                 kNCandSteps };

  struct {
    std::array<int32_t, kNCandSteps> candstats;
    std::array<int32_t, kNCandSteps> truecandstats;
  } statisticsRegistry;

  void resetHistos()
  {
    for (Int_t ii = 0; ii < kNCandSteps; ii++) {
      statisticsRegistry.candstats[ii] = 0;
      statisticsRegistry.truecandstats[ii] = 0;
    }
  }
  void FillCandCounter(int kn, bool istrue = false)
  {
    statisticsRegistry.candstats[kn]++;
    if (istrue) {
      statisticsRegistry.truecandstats[kn]++;
    }
  }
  void fillHistos()
  {
    for (Int_t ii = 0; ii < kNCandSteps; ii++) {
      registry.fill(HIST("hSelectedCandidatesCounter"), ii, statisticsRegistry.candstats[ii]);
      registry.fill(HIST("hSelectedTrueHypertritonCounter"), ii, statisticsRegistry.truecandstats[ii]);
    }
  }

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

    TString CandCounterbinLabel[11] = {"Total", "VtxCosPA", "TrackEta", "MomRapidity", "Lifetime", "VtxDcaDau", "d TOFPID", "TPCPID&Mass", "TPCNcls", "DauPt", "PionDcatoPV"};
    for (int i{0}; i < kNCandSteps; i++) {
      registry.get<TH1>(HIST("hSelectedCandidatesCounter"))->GetXaxis()->SetBinLabel(i + 1, CandCounterbinLabel[i]);
      registry.get<TH1>(HIST("hSelectedTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(i + 1, CandCounterbinLabel[i]);
    }
  }

  //------------------------------------------------------------------
  // Analysis for a single candidate
  template <class TTrackTo, typename TCollisionTable, typename TCandTable>
  void CandidateAnalysis(TCollisionTable const& dCollision, TCandTable const& candData, bool& if_hasvtx, bool isTrueCand = false, double MClifetime = -1, double lPt = -1){

    FillCandCounter(kCandAll, isTrueCand);

      auto track0 = candData.template track0_as<TTrackTo>();
      auto track1 = candData.template track1_as<TTrackTo>();
      auto track2 = candData.template track2_as<TTrackTo>();
 
      if (candData.vtxcosPA(dCollision.posX(), dCollision.posY(), dCollision.posZ()) < vtxcospa) {
        return;
      }
    FillCandCounter(kCandCosPA, isTrueCand);
      if (TMath::Abs(track0.eta()) > etacut || TMath::Abs(track1.eta()) > etacut || TMath::Abs(track2.eta()) > etacut) {
        return;
      }
    FillCandCounter(kCandDauEta, isTrueCand);
      if (TMath::Abs(candData.yHypertriton()) > rapidity) {
        return;
      }
    FillCandCounter(kCandRapidity, isTrueCand);
      double ct = candData.distovertotmom(dCollision.posX(), dCollision.posY(), dCollision.posZ()) * o2::constants::physics::MassHyperTriton;
      if (ct > lifetimecut) {
        return;
      }
    FillCandCounter(kCandct, isTrueCand);
      if (candData.dcaVtxdaughters() > dcavtxdau) {
        return;
      }
    FillCandCounter(kCandDcaDau, isTrueCand);

      registry.fill(HIST("hDeuteronTOFVsPBeforeTOFCut"), track2.sign() * track2.p(), track2.tofNSigmaDe());
      if (isTrueCand) {
        registry.fill(HIST("hDeuteronTOFVsPBeforeTOFCutSig"), track2.sign() * track2.p(), track2.tofNSigmaDe());
      }
      if ((track2.tofNSigmaDe() < TofPidNsigmaMin || track2.tofNSigmaDe() > TofPidNsigmaMax) && track2.p() > minDeuteronPUseTOF) {
        return;
      }
    FillCandCounter(kCandTOFPID, isTrueCand);
      registry.fill(HIST("hDeuteronTOFVsPAtferTOFCut"), track2.sign() * track2.p(), track2.tofNSigmaDe());
      if (isTrueCand) {
        registry.fill(HIST("hDeuteronTOFVsPAtferTOFCutSig"), track2.sign() * track2.p(), track2.tofNSigmaDe());
      }

      // 3sigma region for Dalitz plot
      double lowersignallimit = o2::constants::physics::MassHyperTriton - 3 * mcsigma;
      double uppersignallimit = o2::constants::physics::MassHyperTriton + 3 * mcsigma;

      // Hypertriton
      if (TMath::Abs(track0.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(track1.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(track2.tpcNSigmaDe()) < TpcPidNsigmaCut && candData.mHypertriton() > h3LMassLowerlimit && candData.mHypertriton() < h3LMassUpperlimit) {

    FillCandCounter(kCandTPCPID, isTrueCand);

        if (track0.tpcNClsCrossedRows() > mincrossedrowsproton && track1.tpcNClsCrossedRows() > mincrossedrowspion && track2.tpcNClsCrossedRows() > mincrossedrowsdeuteron) {

    FillCandCounter(kCandTPCNcls, isTrueCand);

          if (candData.track0pt() > minProtonPt && candData.track0pt() < maxProtonPt && candData.track1pt() > minPionPt && candData.track1pt() < maxPionPt && candData.track2pt() > minDeuteronPt && candData.track2pt() < maxDeuteronPt) {
    FillCandCounter(kCandDauPt, isTrueCand);

            if (TMath::Abs(candData.dcatrack1topv()) > dcapiontopv) {
              if_hasvtx = true;
    FillCandCounter(kCandDcaToPV, isTrueCand);

              registry.fill(HIST("hVtxCosPA"), candData.vtxcosPA(dCollision.posX(), dCollision.posY(), dCollision.posZ()));
              registry.fill(HIST("hDCAVtxDau"), candData.dcaVtxdaughters());
              registry.fill(HIST("hPtProton"), candData.track0pt());
              registry.fill(HIST("hPtAntiPion"), candData.track1pt());
              registry.fill(HIST("hPtDeuteron"), candData.track2pt());
              registry.fill(HIST("hDCAProtonToPV"), candData.dcatrack0topv());
              registry.fill(HIST("hDCAPionToPV"), candData.dcatrack1topv());
              registry.fill(HIST("hDCADeuteronToPV"), candData.dcatrack2topv());
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
              registry.fill(HIST("hProtonTPCVsPt"), candData.track0pt(), track0.tpcNSigmaPr());
              registry.fill(HIST("hPionTPCVsPt"), candData.track1pt(), track1.tpcNSigmaPi());
              registry.fill(HIST("hDeuteronTPCVsPt"), candData.track2pt(), track2.tpcNSigmaDe());
              registry.fill(HIST("hMassHypertriton"), candData.mHypertriton());
              registry.fill(HIST("hMassHypertritonTotal"), candData.mHypertriton());
              registry.fill(HIST("h3dMassHypertriton"), 0., candData.pt(), candData.mHypertriton()); // dCollision.centV0M() instead of 0. once available
              registry.fill(HIST("h3dTotalHypertriton"), ct, candData.pt(), candData.mHypertriton());
              if (candData.mHypertriton() > lowersignallimit && candData.mHypertriton() < uppersignallimit) {
                registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}),
                              RecoDecay::m2(array{array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
              }
              if (isTrueCand) {
                registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, candData.mHypertriton());
              }

              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassHypertritonDca"), candData.dcaVtxdaughters(), candData.pt(), candData.mHypertriton());
              }
            }
          }
        }
      }

      // AntiHypertriton
      if (TMath::Abs(track0.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(track1.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(track2.tpcNSigmaDe()) < TpcPidNsigmaCut && candData.mAntiHypertriton() > h3LMassLowerlimit && candData.mAntiHypertriton() < h3LMassUpperlimit) {

    FillCandCounter(kCandTOFPID, isTrueCand);

        if (track0.tpcNClsCrossedRows() > mincrossedrowspion && track1.tpcNClsCrossedRows() > mincrossedrowsproton && track2.tpcNClsCrossedRows() > mincrossedrowsdeuteron) {

    FillCandCounter(kCandTPCNcls, isTrueCand);

          if (candData.track0pt() > minPionPt && candData.track0pt() < maxPionPt && candData.track1pt() > minProtonPt && candData.track1pt() < maxProtonPt && candData.track2pt() > minDeuteronPt && candData.track2pt() < maxDeuteronPt) {
    FillCandCounter(kCandDauPt, isTrueCand);
            if (TMath::Abs(candData.dcatrack0topv()) > dcapiontopv) {
              if_hasvtx = true;
    FillCandCounter(kCandDcaToPV, isTrueCand);

              registry.fill(HIST("hVtxCosPA"), candData.vtxcosPA(dCollision.posX(), dCollision.posY(), dCollision.posZ()));
              registry.fill(HIST("hDCAVtxDau"), candData.dcaVtxdaughters());
              registry.fill(HIST("hPtAntiProton"), candData.track1pt());
              registry.fill(HIST("hPtPion"), candData.track0pt());
              registry.fill(HIST("hPtAntiDeuteron"), candData.track2pt());
              registry.fill(HIST("hDCAProtonToPV"), candData.dcatrack1topv());
              registry.fill(HIST("hDCAPionToPV"), candData.dcatrack0topv());
              registry.fill(HIST("hDCADeuteronToPV"), candData.dcatrack2topv());
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
              registry.fill(HIST("hProtonTPCVsPt"), candData.track1pt(), track1.tpcNSigmaPr());
              registry.fill(HIST("hPionTPCVsPt"), candData.track0pt(), track0.tpcNSigmaPi());
              registry.fill(HIST("hDeuteronTPCVsPt"), candData.track2pt(), track2.tpcNSigmaDe());
              registry.fill(HIST("hMassAntiHypertriton"), candData.mAntiHypertriton());
              registry.fill(HIST("hMassHypertritonTotal"), candData.mAntiHypertriton());
              registry.fill(HIST("h3dMassAntiHypertriton"), 0., candData.pt(), candData.mAntiHypertriton()); // dCollision.centV0M() instead of 0. once available
              registry.fill(HIST("h3dTotalHypertriton"), ct, candData.pt(), candData.mAntiHypertriton());
              if (candData.mAntiHypertriton() > lowersignallimit && candData.mAntiHypertriton() < uppersignallimit) {
                registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}),
                              RecoDecay::m2(array{array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
              }

              if (isTrueCand) {
                registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, candData.mHypertriton());
              }
              if (saveDcaHist == 1) {
                registry.fill(HIST("h3dMassAntiHypertritonDca"), candData.dcaVtxdaughters(), candData.pt(), candData.mAntiHypertriton());
              }
            }
          }
        }
      }

  }

  //------------------------------------------------------------------
  // process real data analysis
  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, FullTracksExtIU const& tracks)
  {
    registry.fill(HIST("hSelectedEventCounter"), 0.5);
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hSelectedEventCounter"), 1.5);

    bool if_hasvtx = false;

    for (auto& vtx : vtx3bodydatas) {
      CandidateAnalysis<FullTracksExtIU>(collision, vtx, if_hasvtx);
    }

    if (if_hasvtx)
      registry.fill(HIST("hSelectedEventCounter"), 2.5);
    fillHistos();
    resetHistos();
  }
  PROCESS_SWITCH(hypertriton3bodyAnalysis, processData, "Real data analysis", true);

  //------------------------------------------------------------------
  // process mc analysis
  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, aod::McParticles const& mcparticles, MCLabeledFullTracksExtIU const& tracks)
  {
    registry.fill(HIST("hSelectedEventCounter"), 0.5);
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hSelectedEventCounter"), 1.5);

    bool if_hasvtx = false;

    for (auto& vtx : vtx3bodydatas) {
      // int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      double MClifetime = -1;
      bool isTrueCand = false;
      auto track0 = vtx.track0_as<MCLabeledFullTracksExtIU>();
      auto track1 = vtx.track1_as<MCLabeledFullTracksExtIU>();
      auto track2 = vtx.track2_as<MCLabeledFullTracksExtIU>();
      if (track0.has_mcParticle() && track1.has_mcParticle() && track2.has_mcParticle()) {
        auto lMCTrack0 = track0.mcParticle_as<aod::McParticles>();
        auto lMCTrack1 = track1.mcParticle_as<aod::McParticles>();
        auto lMCTrack2 = track2.mcParticle_as<aod::McParticles>();
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
                    isTrueCand = true;
                    MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
                  }
                }
              }
            }
          }
        }
      }

      CandidateAnalysis<MCLabeledFullTracksExtIU>(collision, vtx, if_hasvtx, isTrueCand, MClifetime, lPt);
    }

    if (if_hasvtx)
      registry.fill(HIST("hSelectedEventCounter"), 2.5);
    fillHistos();
    resetHistos();
  }
  PROCESS_SWITCH(hypertriton3bodyAnalysis, processMC, "MC analysis", false);

};

// check vtx3body with mclabels
struct hypertriton3bodyLabelCheck {
  HistogramRegistry registry{
    "registry",
    {
      {"hLabeledVtxCounter", "hLabeledVtxCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hMassTrueH3L", "hMassTrueH3L", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
      {"hMassTrueH3LMatter", "hMassTrueH3LMatter", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
      {"hMassTrueH3LAntiMatter", "hMassTrueH3LAntiMatter", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
    },
  };

  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("hLabeledVtxCounter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hLabeledVtxCounter"))->GetXaxis()->SetBinLabel(2, "TrueMCH3L");
    registry.get<TH1>(HIST("hLabeledVtxCounter"))->GetXaxis()->SetBinLabel(3, "Nonrepetitive");
  }

  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision)
  {
    // dummy function
  }
  PROCESS_SWITCH(hypertriton3bodyLabelCheck, process, "Donot check MC label tables", true);

  void processCheckLabel(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::Vtx3BodyDatas, aod::McVtx3BodyLabels> const& vtx3bodydatas, MCLabeledFullTracksExtIU const& tracks, aod::McParticles const& particlesMC)
  {
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }

    std::vector<int64_t> set_mothertrack;
    for (auto& vtx : vtx3bodydatas) {
      registry.fill(HIST("hLabeledVtxCounter"), 0.5);
      if (vtx.mcParticleId() != -1) {
        auto mcparticle = vtx.mcParticle_as<aod::McParticles>();
        if (mcparticle.pdgCode() == 1010010030) {
          registry.fill(HIST("hLabeledVtxCounter"), 1.5);
          registry.fill(HIST("hMassTrueH3L"), vtx.mHypertriton());
          registry.fill(HIST("hMassTrueH3LMatter"), vtx.mHypertriton());
          auto p = std::find(set_mothertrack.begin(), set_mothertrack.end(), mcparticle.globalIndex());
          if (p == set_mothertrack.end()) {
            set_mothertrack.push_back(mcparticle.globalIndex());
            registry.fill(HIST("hLabeledVtxCounter"), 2.5);
          }
        } else if (mcparticle.pdgCode() == -1010010030) {
          registry.fill(HIST("hLabeledVtxCounter"), 1.5);
          registry.fill(HIST("hMassTrueH3L"), vtx.mAntiHypertriton());
          registry.fill(HIST("hMassTrueH3LAntiMatter"), vtx.mAntiHypertriton());
          auto p = std::find(set_mothertrack.begin(), set_mothertrack.end(), mcparticle.globalIndex());
          if (p == set_mothertrack.end()) {
            set_mothertrack.push_back(mcparticle.globalIndex());
            registry.fill(HIST("hLabeledVtxCounter"), 2.5);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(hypertriton3bodyLabelCheck, processCheckLabel, "Check MC label tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertriton3bodyAnalysis>(cfgc),
    adaptAnalysisTask<hypertriton3bodyQa>(cfgc),
    adaptAnalysisTask<hypertriton3bodyLabelCheck>(cfgc),
  };
}
