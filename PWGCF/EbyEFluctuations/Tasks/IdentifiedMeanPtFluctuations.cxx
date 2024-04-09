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

/// \author Sweta Singh (sweta.singh@cern.ch)

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

namespace o2::aod
{
using MyCollisions = soa::Join<aod::Collisions,
                               aod::EvSels,
                               aod::Mults,
                               aod::CentFT0Cs>;
using MyTracks = soa::Join<aod::FullTracks,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                           aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::StoredTracks,
                           aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFbeta, aod::TOFSignal,
                           aod::TracksExtra, aod::TracksIU, aod::TracksDCA, aod::TrackSelection>;

using MyCollision = MyCollisions::iterator;
using MyTrack = MyTracks::iterator;
} // namespace o2::aod

double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
double massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
double massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

struct IdentifiedMeanPtFluctuations {

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec vtxZAxis = {100, -20, 20, "Z (cm)"};
    AxisSpec dcaAxis = {1002, -5.01, 5.01, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {1002, -5.01, 5.01, "DCA_{z} (cm)"};
    AxisSpec ptAxis = {400, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis = {400, 0.0, 4.0, "#it{p} (GeV/#it{c})"};
    AxisSpec betaAxis = {200, 0.0, 2.0, "TOF_{#beta} (GeV/#it{c})"};
    AxisSpec dEdxAxis = {2000, 0.0, 200.0, "dE/dx (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -1.5, 1.5, "#eta"};
    AxisSpec nSigmaTPCAxis = {100, -5., 5., "n#sigma_{TPC}^{proton}"};
    AxisSpec nSigmaTPCAxispid = {110, -5.5, 5.5, "n#sigma_{TPC}"};
    AxisSpec nSigmaTOFAxispid = {110, -5.5, 5.5, "n#sigma_{TOF}"};
    // AxisSpec nChAxis = {2500, -0.5, 2499.5, "nCh"};
    AxisSpec centAxis = {100, 0., 100., "centrality"};
    AxisSpec subAxis = {30, 0., 30., "sample"};
    AxisSpec nchAxis = {3200, 0., 3200., "nch"};
    AxisSpec varAxis1 = {400, 0., 4., "var1"};
    AxisSpec varAxis2 = {400, 0., 4., "var2"};
    AxisSpec Chi2Axis = {100, 0., 100., "Chi2"};
    AxisSpec CrossedrowTPCAxis = {600, 0., 600., "TPC Crossed rows"};
    AxisSpec Counter = {10, 0., 10., "events"};

    // QA Plots
    histos.add("hEventCounter", "event counts", kTH1D, {Counter});

    auto h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1D, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Global track passed");
    h->GetXaxis()->SetBinLabel(3, "DCAxy passed");
    h->GetXaxis()->SetBinLabel(4, "DCAz passed");
    h->GetXaxis()->SetBinLabel(5, "Eta-cut passed");
    h->GetXaxis()->SetBinLabel(6, "pT-cut passed");
    h->GetXaxis()->SetBinLabel(7, "TPC crossed rows passed");
    h->GetXaxis()->SetBinLabel(8, "TPC Chai2cluster passed");
    h->GetXaxis()->SetBinLabel(9, "ITS Chai2cluster passed");

    histos.add("hZvtx_before_sel", "hZvtx_before_sel", kTH1D, {vtxZAxis});
    histos.add("hZvtx_after_sel", "hZvtx_after_sel", kTH1D, {vtxZAxis});
    histos.add("hZvtx_after_sel8", "hZvtx_after_sel8", kTH1D, {vtxZAxis});
    histos.add("hP", "hP", kTH1D, {pAxis});
    histos.add("hEta", ";hEta", kTH1D, {etaAxis});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("hNsigmaTPC", "hNsigmaTPC", kTH2D,
               {pAxis, nSigmaTPCAxis});
    histos.add("hDCAxy", "hDCAxy", kTH1D, {dcaAxis});
    histos.add("hDCAz", "hDCAz", kTH1D, {dcazAxis});

    histos.add("hPtDCAxy", "hPtDCAxy", kTH2D, {ptAxis, dcaAxis});
    histos.add("hPtDCAz", "hPtDCAz", kTH2D, {ptAxis, dcazAxis});
    histos.add("NSigamaTPCpion", "NSigamaTPCpion", kTH2D, {ptAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCkaon", "NSigamaTPCkaon", kTH2D, {ptAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCproton", "NSigamaTPCproton", kTH2D, {ptAxis, nSigmaTPCAxispid});

    histos.add("NSigamaTOFpion", "NSigamaTOFpion", kTH2D, {ptAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFkaon", "NSigamaTOFkaon", kTH2D, {ptAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFproton", "NSigamaTOFproton", kTH2D, {ptAxis, nSigmaTOFAxispid});

    histos.add("NSigamaTPCTOFpion", "NSigamaTPCTOFpion", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFkaon", "NSigamaTPCTOFkaon", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFproton", "NSigamaTPCTOFproton", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});

    histos.add("hPtPion", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("hPtKaon", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("hPtProton", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});

    histos.add("hEtaPion", ";hEta", kTH1D, {etaAxis});
    histos.add("hEtaKaon", ";hEta", kTH1D, {etaAxis});
    histos.add("hEtaProton", ";hEta", kTH1D, {etaAxis});
    //=====================rapidity=====================================
    histos.add("hyPion", ";hyPion", kTH1D, {etaAxis});
    histos.add("hyKaon", ";hyKaon", kTH1D, {etaAxis});
    histos.add("hyProton", ";hyProton", kTH1D, {etaAxis});

    histos.add("hPtCh", "hPtCh", kTH2D, {nchAxis, ptAxis});
    histos.add("hPtChPion", "hPtChPion", kTH2D, {nchAxis, ptAxis});
    histos.add("hPtChKaon", "hPtChKaon", kTH2D, {nchAxis, ptAxis});
    histos.add("hPtChProton", "hPtChProton", kTH2D, {nchAxis, ptAxis});

    histos.add("hPtCent", "hPtCent", kTH2D, {centAxis, ptAxis});
    histos.add("hPtCentPion", "hPtCentPion", kTH2D, {centAxis, ptAxis});
    histos.add("hPtCentKaon", "hPtCentKaon", kTH2D, {centAxis, ptAxis});
    histos.add("hPtCentProton", "hPtCentProton", kTH2D, {centAxis, ptAxis});

    histos.add("hMeanPtCh", "hMeanPtCh", kTH2D, {nchAxis, ptAxis});
    histos.add("hCent", "hCent", kTH2D, {nchAxis, centAxis});

    histos.add("hVar1", "hVar1", kTH2D, {subAxis, centAxis});
    histos.add("hVar2", "hVar2", kTH2D, {subAxis, centAxis});
    histos.add("hVar2meanpt", "hVar2meanpt", kTH2D, {centAxis, varAxis2});
    histos.add("hVar", "hVar", kTH2D, {subAxis, centAxis});
    histos.add("hVarc", "hVarc", kTH2D, {subAxis, centAxis});

    histos.add("hVar1pi", "hVar1pi", kTH2D, {subAxis, centAxis});
    histos.add("hVar2pi", "hVar2pi", kTH2D, {subAxis, centAxis});
    histos.add("hVarpi", "hVarpi", kTH2D, {subAxis, centAxis});
    histos.add("hVar2meanptpi", "hVar2meanptpi", kTH2D, {centAxis, varAxis2});

    histos.add("hVar1k", "hVar1k", kTH2D, {subAxis, centAxis});
    histos.add("hVar2k", "hVar2k", kTH2D, {subAxis, centAxis});
    histos.add("hVark", "hVark", kTH2D, {subAxis, centAxis});
    histos.add("hVar2meanptk", "hVar2meanptk", kTH2D, {centAxis, varAxis2});

    histos.add("hVar1p", "hVar1p", kTH2D, {subAxis, centAxis});
    histos.add("hVar2p", "hVar2p", kTH2D, {subAxis, centAxis});
    histos.add("hVarp", "hVarp", kTH2D, {subAxis, centAxis});
    histos.add("hVar2meanptp", "hVar2meanptp", kTH2D, {centAxis, varAxis2});

    //--------------------------------nch----------------------------------
    histos.add("hVar1x", "hVar1x", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2x", "hVar2x", kTH2D, {subAxis, nchAxis});
    histos.add("hVarx", "hVarx", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptx", "hVar2meanptx", kTH2D, {nchAxis, varAxis2});

    histos.add("hVar1pix", "hVar1pix", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2pix", "hVar2pix", kTH2D, {subAxis, nchAxis});
    histos.add("hVarpix", "hVarpix", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptpix", "hVar2meanptpix", kTH2D, {nchAxis, varAxis2});

    histos.add("hVar1kx", "hVar1kx", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2kx", "hVar2kx", kTH2D, {subAxis, nchAxis});
    histos.add("hVarkx", "hVarkx", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptkx", "hVar2meanptkx", kTH2D, {nchAxis, varAxis2});

    histos.add("hVar1px", "hVar1px", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2px", "hVar2px", kTH2D, {subAxis, nchAxis});
    histos.add("hVarpx", "hVarpx", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptpx", "hVar2meanptpx", kTH2D, {nchAxis, varAxis2});

    histos.add("ht", "ht", kTH1D, {centAxis});

    histos.add("hCentrality", "hCentrality", kTH1D, {centAxis});

    histos.add("hPEta", "hPEta", kTH2D, {pAxis, etaAxis});
    histos.add("hPtEta", "hPtEta", kTH2D, {ptAxis, etaAxis});
    histos.add("hPy", "hPy", kTH2D, {pAxis, etaAxis});
    histos.add("hPty", "hPty", kTH2D, {ptAxis, etaAxis});

    histos.add("hPtyPion", "hPtyPion", kTH2D, {ptAxis, etaAxis});
    histos.add("hPtyKaon", "hPtyKaon", kTH2D, {ptAxis, etaAxis});
    histos.add("hPtyProton", "hPtyProton", kTH2D, {ptAxis, etaAxis});

    histos.add("hTOFbeta", "hTOFbeta", kTH2D, {pAxis, betaAxis});
    histos.add("hdEdx", "hdEdx", kTH2D, {pAxis, dEdxAxis});

    histos.add("hTOFbeta_afterselection", "hTOFbeta_afterselection", kTH2D, {pAxis, betaAxis});
    histos.add("hdEdx_afterselection", "hdEdx_afterselection", kTH2D, {pAxis, dEdxAxis});

    histos.add("hTOFbeta_afterselection1", "hTOFbeta_afterselection1", kTH2D, {pAxis, betaAxis});
    histos.add("hdEdx_afterselection1", "hdEdx_afterselection1", kTH2D, {pAxis, dEdxAxis});

    histos.add("hTPCchi2perCluster_before", "TPC #Chi^{2}/Cluster", kTH1D, {Chi2Axis});
    histos.add("hITSchi2perCluster_before", "ITS #Chi^{2}/Cluster", kTH1D, {Chi2Axis});
    histos.add("hTPCCrossedrows_before", "Crossed TPC rows", kTH1D, {CrossedrowTPCAxis});

    histos.add("hTPCchi2perCluster_after", "TPC #Chi^{2}/Cluster", kTH1D, {Chi2Axis});
    histos.add("hITSchi2perCluster_after", "ITS #Chi^{2}/Cluster", kTH1D, {Chi2Axis});
    histos.add("hTPCCrossedrows_after", "Crossed TPC rows", kTH1D, {CrossedrowTPCAxis});
  }

  void process(aod::MyCollision const& coll, aod::MyTracks const& inputTracks)

  {
    histos.fill(HIST("hEventCounter"), 1.);

    histos.fill(HIST("hZvtx_before_sel"), coll.posZ());
    if (fabs(coll.posZ()) > 10.f) {
      return;
    }

    histos.fill(HIST("hEventCounter"), 2.);

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());

    if (!coll.sel8()) {
      return;
    }
    histos.fill(HIST("hZvtx_after_sel8"), coll.posZ());

    histos.fill(HIST("hEventCounter"), 3.);

    const auto cent = coll.centFT0C();
    histos.fill(HIST("hCentrality"), cent);

    double nCh = 0.;
    double nChpi = 0.;
    double nChk = 0.;
    double nChp = 0.;

    double Q1 = 0, Q2 = 0;
    double Q1pi = 0, Q2pi = 0;
    double Q1k = 0, Q2k = 0;
    double Q1p = 0, Q2p = 0;
    double var1 = 0, var2 = 0, twopar_allcharge = 0;
    double var1pi = 0, var2pi = 0;
    double var1k = 0, var2k = 0;
    double var1p = 0, var2p = 0;
    //  cent = 0;

    // sampling
    int sample = histos.get<TH1>(HIST("hZvtx_after_sel8"))->GetEntries();
    sample = sample % 15;

    // Perfroming the track  selection==========================================
    for (auto track : inputTracks) {
      // Loop over tracks

      // inital tracks
      histos.fill(HIST("tracksel"), 1);

      histos.fill(HIST("hTPCchi2perCluster_before"), track.tpcChi2NCl());
      histos.fill(HIST("hITSchi2perCluster_before"), track.itsChi2NCl());
      histos.fill(HIST("hTPCCrossedrows_before"), track.tpcNClsCrossedRows());

      // tracks passed after GlobalTrackcut
      if (!track.isGlobalTrack())
        continue;
      histos.fill(HIST("tracksel"), 2);

      // tracks passed after DCAxy
      //        if (!(fabs(track.dcaXY()) < 0.12)) continue;//global cut already includes
      histos.fill(HIST("tracksel"), 3);

      // tracks passed after DCAz
      //
      histos.fill(HIST("hDCAxy"), track.dcaXY());
      histos.fill(HIST("hDCAz"), track.dcaZ());

      //    if (!(fabs(track.dcaZ()) < 1.)) continue;//global cut already includes (DCAz< 2.0) cm
      histos.fill(HIST("tracksel"), 4);

      // tracks passed after Eta-cut
      if (!(fabs(track.eta()) < 0.8))
        continue;
      histos.fill(HIST("tracksel"), 5);

      // tracks passed after pT-cut
      if (!(track.pt() > 0.0 && track.pt() < 2.))
        continue; // pt = 0.15
      histos.fill(HIST("tracksel"), 6);

      //    if (track.tpcNClsCrossedRows() < 70.0) continue;
      histos.fill(HIST("hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
      histos.fill(HIST("tracksel"), 7);

      //      if (track.tpcChi2NCl() > 4.0) continue;
      histos.fill(HIST("hTPCchi2perCluster_after"), track.tpcChi2NCl());
      histos.fill(HIST("tracksel"), 8);

      //      if (track.itsChi2NCl() > 36.0) continue;
      histos.fill(HIST("hITSchi2perCluster_after"), track.itsChi2NCl());
      histos.fill(HIST("tracksel"), 9);

      nCh += 1.;

      Q1 += track.pt();
      Q2 += (track.pt() * track.pt());

      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPtDCAxy"), track.pt(), track.dcaXY());
      histos.fill(HIST("hPtDCAz"), track.pt(), track.dcaZ());

      histos.fill(HIST("hPtEta"), track.pt(), track.eta());
      histos.fill(HIST("hPEta"), track.p(), track.eta());

      histos.fill(HIST("hNsigmaTPC"), track.p(), track.tpcNSigmaPr());

      // only TPC tracks: Pion, Kaon, Proton
      if (track.hasTPC() && abs(track.tpcNSigmaPi()) < 2.)
        histos.fill(HIST("NSigamaTPCpion"), track.pt(), track.tpcNSigmaPi());
      if (track.hasTPC() && abs(track.tpcNSigmaKa()) < 2.)
        histos.fill(HIST("NSigamaTPCkaon"), track.pt(), track.tpcNSigmaKa());
      if (track.hasTPC() && abs(track.tpcNSigmaPr()) < 2.)
        histos.fill(HIST("NSigamaTPCproton"), track.pt(), track.tpcNSigmaPr());

      // only TOF tracks: Pion, Kaon, Proton
      if (track.hasTOF() && abs(track.tofNSigmaPi()) < 2.)
        histos.fill(HIST("NSigamaTOFpion"), track.pt(), track.tofNSigmaPi());
      if (track.hasTOF() && abs(track.tofNSigmaKa()) < 2.)
        histos.fill(HIST("NSigamaTOFkaon"), track.pt(), track.tofNSigmaKa());
      if (track.hasTOF() && abs(track.tofNSigmaPr()) < 2.)
        histos.fill(HIST("NSigamaTOFproton"), track.pt(), track.tofNSigmaPr());

      if (track.hasTPC())
        histos.fill(HIST("hdEdx"), track.p(), track.tpcSignal());
      if (track.hasTOF())
        histos.fill(HIST("hTOFbeta"), track.p(), track.beta());

      //=============================pion==============================================================
      // only TPC+TOF tracks: Pion, Kaon, Proton
      if ((track.hasTPC() && abs(track.tpcNSigmaPi()) < 2.) && (track.hasTOF() && abs(track.tofNSigmaPi()) < 2.)) {
        histos.fill(HIST("NSigamaTPCTOFpion"), track.tpcNSigmaPi(), track.tofNSigmaPi());

        histos.fill(HIST("hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection"), track.p(), track.beta());
      }

      // pion-TPC-----------------------------------------------------------------------------------

      if ((track.hasTPC() && abs(track.tpcNSigmaPi()) < 2. && (track.pt() >= 0.15 && track.pt() < 0.65) && (abs(track.rapidity(massPi)) < 0.5) && (std::abs(track.tofNSigmaEl()) > 1.0 && std::abs(track.tpcNSigmaKa()) > 2.0 && std::abs(track.tpcNSigmaPr()) > 2.0))) {

        histos.fill(HIST("hPtPion"), track.pt());
        histos.fill(HIST("hEtaPion"), track.eta());
        histos.fill(HIST("hyPion"), track.rapidity(massPi));
        histos.fill(HIST("hPtyPion"), track.pt(), track.rapidity(massPi));

        nChpi += 1.;
        Q1pi += track.pt();
        Q2pi += (track.pt() * track.pt());

        histos.fill(HIST("hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      // pion->(TPC+TOF)------------------------------------------------------------------------------------
      if ((track.pt() >= 0.65 && track.pt() < 2.0) && (abs(track.rapidity(massPi)) < 0.5) && track.hasTPC() && track.hasTOF() && (std::abs(track.tofNSigmaEl()) > 1.0 && std::abs(track.tofNSigmaKa()) > 2.0 && std::abs(track.tofNSigmaPr()) > 2.0) && abs(sqrt(track.tpcNSigmaPi()) * (track.tpcNSigmaPi()) + (track.tofNSigmaPi()) * (track.tofNSigmaPi())) < 2.) {

        histos.fill(HIST("hPtPion"), track.pt());
        histos.fill(HIST("hEtaPion"), track.eta());
        histos.fill(HIST("hyPion"), track.rapidity(massPi));
        histos.fill(HIST("hPtyPion"), track.pt(), track.rapidity(massPi));

        nChpi += 1.;
        Q1pi += track.pt();
        Q2pi += (track.pt() * track.pt());

        histos.fill(HIST("hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      //===========================kaon===============================================================

      if ((track.hasTPC() && abs(track.tpcNSigmaKa()) < 2.) && (track.hasTOF() && abs(track.tofNSigmaKa()) < 2.)) {
        histos.fill(HIST("NSigamaTPCTOFkaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());
        histos.fill(HIST("hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (track.hasTPC() && abs(track.tpcNSigmaKa()) < 2. && (track.pt() >= 0.15 && track.pt() < 0.65) && (abs(track.rapidity(massKa)) < 0.5) && (std::abs(track.tofNSigmaEl()) > 1.0 && std::abs(track.tpcNSigmaPi()) > 2.0 && std::abs(track.tpcNSigmaPr()) > 2.0)) {

        histos.fill(HIST("hPtKaon"), track.pt());
        histos.fill(HIST("hEtaKaon"), track.eta());
        histos.fill(HIST("hyKaon"), track.rapidity(massKa));
        histos.fill(HIST("hPtyKaon"), track.pt(), track.rapidity(massKa));

        nChk += 1.;
        Q1k += track.pt();
        Q2k += (track.pt() * track.pt());

        histos.fill(HIST("hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      if ((track.pt() >= 0.65 && track.pt() < 2.0) && (abs(track.rapidity(massKa)) < 0.5) && track.hasTPC() && track.hasTOF() && (std::abs(track.tofNSigmaEl()) > 1.0 && std::abs(track.tofNSigmaPi()) > 2.0 && std::abs(track.tofNSigmaPr()) > 2.0) && (abs(sqrt(track.tpcNSigmaKa()) * (track.tpcNSigmaKa()) + (track.tofNSigmaKa()) * (track.tofNSigmaKa())) < 2.)) {

        histos.fill(HIST("hPtKaon"), track.pt());
        histos.fill(HIST("hEtaKaon"), track.eta());
        histos.fill(HIST("hyKaon"), track.rapidity(massKa));
        histos.fill(HIST("hPtyKaon"), track.pt(), track.rapidity(massKa));

        nChk += 1.;
        Q1k += track.pt();
        Q2k += (track.pt() * track.pt());

        histos.fill(HIST("hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      //============================proton===========================================================

      if ((track.hasTPC() && abs(track.tpcNSigmaPr()) < 2.) && (track.hasTOF() && abs(track.tofNSigmaPr()) < 2.)) {
        histos.fill(HIST("NSigamaTPCTOFproton"), track.tpcNSigmaPr(), track.tofNSigmaPr());

        histos.fill(HIST("hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (track.hasTPC() && abs(track.tpcNSigmaPr()) < 2. && (track.pt() >= 0.4 && track.pt() < 0.85) && (abs(track.rapidity(massPr)) < 0.5) && (std::abs(track.tofNSigmaEl()) > 1.0 && std::abs(track.tpcNSigmaKa()) > 2.0 && std::abs(track.tpcNSigmaPi()) > 2.0)) {

        histos.fill(HIST("hPtProton"), track.pt());
        histos.fill(HIST("hEtaProton"), track.eta());
        histos.fill(HIST("hyProton"), track.rapidity(massPr));
        histos.fill(HIST("hPtyProton"), track.pt(), track.rapidity(massPr));

        nChp += 1.;
        Q1p += track.pt();
        Q2p += (track.pt() * track.pt());

        histos.fill(HIST("hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      if ((track.pt() >= 0.85 && track.pt() < 2.0) && (abs(track.rapidity(massPr)) < 0.5) && track.hasTPC() && track.hasTOF() && (std::abs(track.tofNSigmaEl()) > 1.0 && std::abs(track.tofNSigmaKa()) > 2.0 && std::abs(track.tofNSigmaPi()) > 2.0) && (abs(sqrt(track.tpcNSigmaPr()) * (track.tpcNSigmaPr()) + (track.tofNSigmaPr()) * (track.tofNSigmaPr())) < 2.)) {

        histos.fill(HIST("hPtProton"), track.pt());
        histos.fill(HIST("hEtaProton"), track.eta());
        histos.fill(HIST("hyProton"), track.rapidity(massPr));
        histos.fill(HIST("hPtyProton"), track.pt(), track.rapidity(massPr));

        nChp += 1.;
        Q1p += track.pt();
        Q2p += (track.pt() * track.pt());

        histos.fill(HIST("hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      //=================================================================================================
    }
    // Track loop ends!

    if (nCh < 2)
      return;

    //------------------ all charges-------------------------------------
    var1 = (Q1 * Q1 - Q2) / (nCh * (nCh - 1));
    histos.fill(HIST("hVar1"), sample, cent, var1);
    var2 = (Q1 / nCh);
    histos.fill(HIST("hVar2"), sample, cent, var2);
    histos.fill(HIST("hVarc"), sample, cent);
    histos.fill(HIST("hVar2meanpt"), cent, var2);

    twopar_allcharge = (var1 - var2);
    histos.fill(HIST("hVar"), nCh, twopar_allcharge);

    //---------------------- pions ----------------------------------------

    if (nChpi > 2) {
      var1pi = (Q1pi * Q1pi - Q2pi) / (nChpi * (nChpi - 1));
      histos.fill(HIST("hVar1pi"), sample, cent, var1pi);
      var2pi = (Q1pi / nChpi);
      histos.fill(HIST("hVar2pi"), sample, cent, var2pi);
      histos.fill(HIST("hVar2meanptpi"), cent, var2pi);
    }

    //----------------------- kaons ---------------------------------------
    if (nChk > 2) {
      var1k = (Q1k * Q1k - Q2k) / (nChk * (nChk - 1));
      histos.fill(HIST("hVar1k"), sample, cent, var1k);
      var2k = (Q1k / nChk);
      histos.fill(HIST("hVar2k"), sample, cent, var2k);
      histos.fill(HIST("hVar2meanptk"), cent, var2k);
    }

    //---------------------------- protons ----------------------------------
    if (nChp > 2) {
      var1p = (Q1p * Q1p - Q2p) / (nChp * (nChp - 1));
      histos.fill(HIST("hVar1p"), sample, cent, var1p);
      var2p = (Q1p / nChp);
      histos.fill(HIST("hVar2p"), sample, cent, var2p);
      histos.fill(HIST("hVar2meanptp"), cent, var2p);
    }

    //-----------------------nch-------------------------------------
    histos.fill(HIST("hVar1x"), sample, nCh, var1);
    histos.fill(HIST("hVar2x"), sample, nCh, var2);
    histos.fill(HIST("hVarx"), sample, nCh);
    histos.fill(HIST("hVar2meanptx"), nCh, var2);

    histos.fill(HIST("hVar1pix"), sample, nCh, var1pi);
    histos.fill(HIST("hVar2pix"), sample, nCh, var2pi);
    histos.fill(HIST("hVarpix"), sample, nChpi);
    histos.fill(HIST("hVar2meanptpix"), nCh, var2pi);

    histos.fill(HIST("hVar1kx"), sample, nCh, var1k);
    histos.fill(HIST("hVar2kx"), sample, nCh, var2k);
    histos.fill(HIST("hVarkx"), sample, nChk);
    histos.fill(HIST("hVar2meanptkx"), nCh, var2k);

    histos.fill(HIST("hVar1px"), sample, nCh, var1p);
    histos.fill(HIST("hVar2px"), sample, nCh, var2p);
    histos.fill(HIST("hVarpx"), sample, nChp);
    histos.fill(HIST("hVar2meanptpx"), nCh, var2p);

  } // event loop ends!
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<IdentifiedMeanPtFluctuations>(cfgc)};
  return workflow;
}
