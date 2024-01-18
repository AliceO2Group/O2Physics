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

struct IdentifiedMeanPtFluctuations {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec vtxZAxis = {100, -20, 20, "Z (cm)"};
    AxisSpec dcaAxis = {1002, -5.01, 5.01, "DCA_{xy} (cm)"};
    AxisSpec ptAxis = {40, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis = {40, 0.0, 4.0, "#it{p} (GeV/#it{c})"};
    AxisSpec betaAxis = {200, 0.0, 2.0, "TOF_{#beta} (GeV/#it{c})"};
    AxisSpec dEdxAxis = {2000, 0.0, 200.0, "dE/dx (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -1.5, 1.5, "#eta"};
    AxisSpec nSigmaTPCAxis = {100, -5., 5., "n#sigma_{TPC}^{proton}"};
    AxisSpec nSigmaTPCAxispid = {110, -5.5, 5.5, "n#sigma_{TPC}"};
    AxisSpec nSigmaTOFAxispid = {110, -5.5, 5.5, "n#sigma_{TOF}"};
    AxisSpec nChAxis = {2500, -0.5, 2499.5, "nCh"};
    AxisSpec centAxis = {20, 0., 100., "centrality"};
    AxisSpec subAxis = {30, 0., 30., "sample"};
    AxisSpec nchAxis = {3200, 0., 3200., "nch"};
    AxisSpec varAxis1 = {6000, 0., 6., "var1"};
    AxisSpec varAxis2 = {600, 0., 6., "var2"};

    // QA Plots
    histos.add("hZvtx_before_sel", "hZvtx_before_sel", kTH1F, {vtxZAxis});
    histos.add("hZvtx_after_sel", "hZvtx_after_sel", kTH1F, {vtxZAxis});
    histos.add("hZvtx_after_sel8", "hZvtx_after_sel8", kTH1F, {vtxZAxis});
    histos.add("hP", "hP", kTH1F, {pAxis});
    histos.add("hEta", ";hEta", kTH1F, {etaAxis});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hNsigmaTPC", "hNsigmaTPC", kTH2F,
               {pAxis, nSigmaTPCAxis});
    histos.add("hPtDCAxy", "hPtDCAxy", kTH2F, {ptAxis, dcaAxis});
    histos.add("hPtDCAz", "hPtDCAz", kTH2F, {ptAxis, dcaAxis});
    histos.add("NSigamaTPCpion", "NSigamaTPCpion", kTH2F, {ptAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCkaon", "NSigamaTPCkaon", kTH2F, {ptAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCproton", "NSigamaTPCproton", kTH2F, {ptAxis, nSigmaTPCAxispid});

    histos.add("NSigamaTOFpion", "NSigamaTOFpion", kTH2F, {ptAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFkaon", "NSigamaTOFkaon", kTH2F, {ptAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFproton", "NSigamaTOFproton", kTH2F, {ptAxis, nSigmaTOFAxispid});

    histos.add("NSigamaTPCTOFpion", "NSigamaTPCTOFpion", kTH2F, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFkaon", "NSigamaTPCTOFkaon", kTH2F, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFproton", "NSigamaTPCTOFproton", kTH2F, {nSigmaTPCAxispid, nSigmaTOFAxispid});

    histos.add("hPtPion", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPtKaon", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPtProton", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});

    histos.add("hEtaPion", ";hEta", kTH1F, {etaAxis});
    histos.add("hEtaKaon", ";hEta", kTH1F, {etaAxis});
    histos.add("hEtaProton", ";hEta", kTH1F, {etaAxis});

    histos.add("hPtCh", "hPtCh", kTH2F, {nchAxis, ptAxis});
    histos.add("hPtChPion", "hPtChPion", kTH2F, {nchAxis, ptAxis});
    histos.add("hPtChKaon", "hPtChKaon", kTH2F, {nchAxis, ptAxis});
    histos.add("hPtChProton", "hPtChProton", kTH2F, {nchAxis, ptAxis});

    histos.add("hMeanPtCh", "hMeanPtCh", kTH2F, {nChAxis, ptAxis});
    histos.add("hCent", "hCent", kTH2F, {nChAxis, centAxis});

    histos.add("hVar1", "hVar1", kTH2F, {centAxis, varAxis1});
    histos.add("hVar2", "hVar2", kTH2F, {centAxis, varAxis2});
    histos.add("hVar", "hVar", kTH2F, {subAxis, centAxis});

    histos.add("hVar1pi", "hVar1pi", kTH2F, {centAxis, varAxis1});
    histos.add("hVar2pi", "hVar2pi", kTH2F, {centAxis, varAxis2});
    histos.add("hVarpi", "hVarpi", kTH2F, {subAxis, centAxis});

    histos.add("hVar1k", "hVar1k", kTH2F, {centAxis, varAxis1});
    histos.add("hVar2k", "hVar2k", kTH2F, {centAxis, varAxis2});
    histos.add("hVark", "hVark", kTH2F, {subAxis, centAxis});

    histos.add("hVar1p", "hVar1p", kTH2F, {centAxis, varAxis1});
    histos.add("hVar2p", "hVar2p", kTH2F, {centAxis, varAxis2});
    histos.add("hVarp", "hVarp", kTH2F, {subAxis, centAxis});

    //--------------------------------nch----------------------------------
    histos.add("hVar1x", "hVar1x", kTH2F, {nchAxis, varAxis1});
    histos.add("hVar2x", "hVar2x", kTH2F, {nchAxis, varAxis2});
    histos.add("hVarx", "hVarx", kTH2F, {subAxis, nchAxis});

    histos.add("hVar1pix", "hVar1pix", kTH2F, {nchAxis, varAxis1});
    histos.add("hVar2pix", "hVar2pix", kTH2F, {nchAxis, varAxis2});
    histos.add("hVarpix", "hVarpix", kTH2F, {subAxis, nchAxis});

    histos.add("hVar1kx", "hVar1kx", kTH2F, {nchAxis, varAxis1});
    histos.add("hVar2kx", "hVar2kx", kTH2F, {nchAxis, varAxis2});
    histos.add("hVarkx", "hVarkx", kTH2F, {subAxis, nchAxis});

    histos.add("hVar1px", "hVar1px", kTH2F, {nchAxis, varAxis1});
    histos.add("hVar2px", "hVar2px", kTH2F, {nchAxis, varAxis2});
    histos.add("hVarpx", "hVarpx", kTH2F, {subAxis, nchAxis});

    histos.add("ht", "ht", kTH1F, {centAxis});

    histos.add("hCentrality", "hCentrality", kTH1F, {centAxis});

    histos.add("hPEta", "hPEta", kTH2F, {pAxis, etaAxis});
    histos.add("hPtEta", "hPtEta", kTH2F, {ptAxis, etaAxis});
    histos.add("hPy", "hPy", kTH2F, {pAxis, etaAxis});
    histos.add("hPty", "hPty", kTH2F, {ptAxis, etaAxis});

    histos.add("hTOFbeta", "hTOFbeta", kTH2F, {pAxis, betaAxis});
    histos.add("hdEdx", "hdEdx", kTH2F, {pAxis, dEdxAxis});
  }

  void process(aod::MyCollision const& coll, aod::MyTracks const& inputTracks)

  {

    histos.fill(HIST("hZvtx_before_sel"), coll.posZ());
    if (fabs(coll.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());

    if (!coll.sel8()) {
      return;
    }
    histos.fill(HIST("hZvtx_after_sel8"), coll.posZ());

    const auto cent = coll.centFT0C();
    histos.fill(HIST("hCentrality"), cent);

    float nCh = 0.;
    float nChpi = 0.;
    float nChk = 0.;
    float nChp = 0.;
    std::vector<float> VMeanPt;
    std::vector<float> VMeanPtPion;
    std::vector<float> VMeanPtKaon;
    std::vector<float> VMeanPtProton;

    float Q1 = 0, Q2 = 0;
    float Q1pi = 0, Q2pi = 0;
    float Q1k = 0, Q2k = 0;
    float Q1p = 0, Q2p = 0;
    float var1, var2;
    float var1pi, var2pi;
    float var1k, var2k;
    float var1p, var2p;

    int sample = histos.get<TH1>(HIST("hZvtx_before_sel"))->GetEntries();
    sample = sample % 30;

    // Perfroming the track  selection==============================================================================================================
    for (auto track : inputTracks) {
      // Loop over tracks

      if (!track.isGlobalTrack())
        return;

      if (!((fabs(track.eta()) < 0.8) && (fabs(track.dcaXY()) < 0.12) && (fabs(track.dcaZ()) < 1.) && (track.pt() > 0.15 && track.pt() < 2.))) {
        continue;
      }

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

      VMeanPt.push_back(track.pt());

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

      // only TPC+TOF tracks: Pion, Kaon, Proton
      if ((track.hasTPC() && abs(track.tpcNSigmaPi()) < 2.) && (track.hasTOF() && abs(track.tofNSigmaPi()) < 2.))
        histos.fill(HIST("NSigamaTPCTOFpion"), track.tpcNSigmaPi(), track.tofNSigmaPi());

      if ((track.hasTPC() && abs(track.tpcNSigmaPi()) < 2. && (track.pt() > 0.2 && track.pt() < 0.5)))
        histos.fill(HIST("hPtPion"), track.pt());

      {
        if ((track.hasTPC() && (track.pt() >= 0.5 && track.pt() < 2.0)) && (track.hasTOF() && abs(sqrt(track.tpcNSigmaPi()) * (track.tpcNSigmaPi()) + (track.tofNSigmaPi()) * (track.tofNSigmaPi())) < 2.)) {

          histos.fill(HIST("hPtPion"), track.pt());
          histos.fill(HIST("hEtaPion"), track.eta());

          VMeanPtPion.push_back(track.pt());

          for (ULong64_t jPi = 0; jPi < VMeanPtPion.size(); jPi++) {
            histos.fill(HIST("hPtChPion"), nCh, VMeanPtPion[jPi]);
          }

          nChpi += 1.;
          Q1pi += track.pt();
          Q2pi += (track.pt() * track.pt());
        }
      }

      if ((track.hasTPC() && abs(track.tpcNSigmaKa()) < 2.) && (track.hasTOF() && abs(track.tofNSigmaKa()) < 2.))
        histos.fill(HIST("NSigamaTPCTOFkaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());

      if ((track.hasTPC() && abs(track.tpcNSigmaKa()) < 2. && (track.pt() > 0.2 && track.pt() < 0.5)))
        histos.fill(HIST("hPtKaon"), track.pt());
      {
        if ((track.hasTPC() && (track.pt() >= 0.5 && track.pt() < 2.0)) && (track.hasTOF() && abs(sqrt(track.tpcNSigmaKa()) * (track.tpcNSigmaKa()) + (track.tofNSigmaKa()) * (track.tofNSigmaKa())) < 2.)) {

          histos.fill(HIST("hPtKaon"), track.pt());
          histos.fill(HIST("hEtaKaon"), track.eta());

          VMeanPtKaon.push_back(track.pt());

          for (ULong64_t jKa = 0; jKa < VMeanPtKaon.size(); jKa++) {
            histos.fill(HIST("hPtChKaon"), nCh, VMeanPtKaon[jKa]);
          }

          nChk += 1.;
          Q1k += track.pt();
          Q2k += (track.pt() * track.pt());
        }
      }

      if ((track.hasTPC() && abs(track.tpcNSigmaPr()) < 2.) && (track.hasTOF() && abs(track.tofNSigmaPr()) < 2.))
        histos.fill(HIST("NSigamaTPCTOFproton"), track.tpcNSigmaPr(), track.tofNSigmaPr());

      if ((track.hasTPC() && abs(track.tpcNSigmaPr()) < 2. && (track.pt() > 0.4 && track.pt() < 0.6)))
        histos.fill(HIST("hPtProton"), track.pt());
      {
        if ((track.hasTPC() && (track.pt() >= 0.6 && track.pt() < 2.0)) && (track.hasTOF() && abs(sqrt(track.tpcNSigmaPr()) * (track.tpcNSigmaPr()) + (track.tofNSigmaPr()) * (track.tofNSigmaPr())) < 2.)) {

          histos.fill(HIST("hPtProton"), track.pt());
          histos.fill(HIST("hEtaProton"), track.eta());

          VMeanPtProton.push_back(track.pt());

          for (ULong64_t jPr = 0; jPr < VMeanPtProton.size(); jPr++) {
            histos.fill(HIST("hPtChProton"), nCh, VMeanPtProton[jPr]);
          }

          nChp += 1.;
          Q1p += track.pt();
          Q2p += (track.pt() * track.pt());
        }
      }
    }
    // Track loop ends!
    VMeanPtPion.clear();
    VMeanPtKaon.clear();
    VMeanPtProton.clear();

    if (nCh < 2)
      return;

    //------------------ all charges-------------------------------------
    var1 = (Q1 * Q1 - Q2) / (nCh * (nCh - 1));
    histos.fill(HIST("hVar1"), cent, var1);
    var2 = (Q1 / nCh);
    histos.fill(HIST("hVar2"), cent, var2);
    histos.fill(HIST("hVar"), sample, cent);

    //---------------------- pions ----------------------------------------

    var1pi = (Q1pi * Q1pi - Q2pi) / (nChpi * (nChpi - 1));
    histos.fill(HIST("hVar1pi"), cent, var1pi);
    var2pi = (Q1pi / nChpi);
    histos.fill(HIST("hVar2pi"), cent, var2pi);

    //----------------------- kaons ---------------------------------------

    var1k = (Q1k * Q1k - Q2k) / (nChk * (nChk - 1));
    histos.fill(HIST("hVar1k"), cent, var1k);
    var2k = (Q1k / nChk);
    histos.fill(HIST("hVar2k"), cent, var2k);

    //---------------------------- protons ----------------------------------

    var1p = (Q1p * Q1p - Q2p) / (nChp * (nChp - 1));
    histos.fill(HIST("hVar1p"), cent, var1p);
    var2p = (Q1p / nChp);
    histos.fill(HIST("hVar2p"), cent, var2p);

    //-----------------------nch-------------------------------------
    histos.fill(HIST("hVar1x"), nCh, var1);
    histos.fill(HIST("hVar2x"), nCh, var2);
    histos.fill(HIST("hVarx"), sample, nCh);

    histos.fill(HIST("hVar1pix"), nCh, var1pi);
    histos.fill(HIST("hVar2pix"), nCh, var2pi);
    histos.fill(HIST("hVarpix"), sample, nChpi);

    histos.fill(HIST("hVar1kx"), nCh, var1k);
    histos.fill(HIST("hVar2kx"), nCh, var2k);
    histos.fill(HIST("hVarkx"), sample, nChk);

    histos.fill(HIST("hVar1px"), nCh, var1p);
    histos.fill(HIST("hVar2px"), nCh, var2p);
    histos.fill(HIST("hVarpx"), sample, nChp);

    for (ULong64_t j = 0; j < VMeanPt.size(); j++) {
      histos.fill(HIST("hPtCh"), nCh, VMeanPt[j]);
    }

    VMeanPt.clear();

  } // event loop ends!
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<IdentifiedMeanPtFluctuations>(cfgc)};
  return workflow;
}
