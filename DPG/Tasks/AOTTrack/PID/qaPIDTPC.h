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

///
/// \file   qaPIDTPC.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Header file for QA tasks of the TPC PID quantities
///

#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TPC QA plots
struct tpcPidQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hdelta[Np] = {"delta/El", "delta/Mu", "delta/Pi",
                                                  "delta/Ka", "delta/Pr", "delta/De",
                                                  "delta/Tr", "delta/He", "delta/Al"};
  static constexpr std::string_view hdelta_pt_pos[Np] = {"delta/pt/pos/El", "delta/pt/pos/Mu", "delta/pt/pos/Pi",
                                                         "delta/pt/pos/Ka", "delta/pt/pos/Pr", "delta/pt/pos/De",
                                                         "delta/pt/pos/Tr", "delta/pt/pos/He", "delta/pt/pos/Al"};
  static constexpr std::string_view hdelta_pt_neg[Np] = {"delta/pt/neg/El", "delta/pt/neg/Mu", "delta/pt/neg/Pi",
                                                         "delta/pt/neg/Ka", "delta/pt/neg/Pr", "delta/pt/neg/De",
                                                         "delta/pt/neg/Tr", "delta/pt/neg/He", "delta/pt/neg/Al"};
  static constexpr std::string_view hexpsigma[Np] = {"expsigma/El", "expsigma/Mu", "expsigma/Pi",
                                                     "expsigma/Ka", "expsigma/Pr", "expsigma/De",
                                                     "expsigma/Tr", "expsigma/He", "expsigma/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigma_pt[Np] = {"nsigma/pt/El", "nsigma/pt/Mu", "nsigma/pt/Pi",
                                                      "nsigma/pt/Ka", "nsigma/pt/Pr", "nsigma/pt/De",
                                                      "nsigma/pt/Tr", "nsigma/pt/He", "nsigma/pt/Al"};
  static constexpr std::string_view hnsigma_pt_pos[Np] = {"nsigma/pt/pos/El", "nsigma/pt/pos/Mu", "nsigma/pt/pos/Pi",
                                                          "nsigma/pt/pos/Ka", "nsigma/pt/pos/Pr", "nsigma/pt/pos/De",
                                                          "nsigma/pt/pos/Tr", "nsigma/pt/pos/He", "nsigma/pt/pos/Al"};
  static constexpr std::string_view hnsigma_pt_neg[Np] = {"nsigma/pt/neg/El", "nsigma/pt/neg/Mu", "nsigma/pt/neg/Pi",
                                                          "nsigma/pt/neg/Ka", "nsigma/pt/neg/Pr", "nsigma/pt/neg/De",
                                                          "nsigma/pt/neg/Tr", "nsigma/pt/neg/He", "nsigma/pt/neg/Al"};

  // With TOF
  static constexpr std::string_view hexpected_wTOF[Np] = {"wTOF/expected/El", "wTOF/expected/Mu", "wTOF/expected/Pi",
                                                          "wTOF/expected/Ka", "wTOF/expected/Pr", "wTOF/expected/De",
                                                          "wTOF/expected/Tr", "wTOF/expected/He", "wTOF/expected/Al"};
  static constexpr std::string_view hdelta_wTOF[Np] = {"wTOF/delta/El", "wTOF/delta/Mu", "wTOF/delta/Pi",
                                                       "wTOF/delta/Ka", "wTOF/delta/Pr", "wTOF/delta/De",
                                                       "wTOF/delta/Tr", "wTOF/delta/He", "wTOF/delta/Al"};
  static constexpr std::string_view hdelta_pt_pos_wTOF[Np] = {"wTOF/delta/pt/pos/El", "wTOF/delta/pt/pos/Mu", "wTOF/delta/pt/pos/Pi",
                                                              "wTOF/delta/pt/pos/Ka", "wTOF/delta/pt/pos/Pr", "wTOF/delta/pt/pos/De",
                                                              "wTOF/delta/pt/pos/Tr", "wTOF/delta/pt/pos/He", "wTOF/delta/pt/pos/Al"};
  static constexpr std::string_view hdelta_pt_neg_wTOF[Np] = {"wTOF/delta/pt/neg/El", "wTOF/delta/pt/neg/Mu", "wTOF/delta/pt/neg/Pi",
                                                              "wTOF/delta/pt/neg/Ka", "wTOF/delta/pt/neg/Pr", "wTOF/delta/pt/neg/De",
                                                              "wTOF/delta/pt/neg/Tr", "wTOF/delta/pt/neg/He", "wTOF/delta/pt/neg/Al"};
  static constexpr std::string_view hexpsigma_wTOF[Np] = {"wTOF/expsigma/El", "wTOF/expsigma/Mu", "wTOF/expsigma/Pi",
                                                          "wTOF/expsigma/Ka", "wTOF/expsigma/Pr", "wTOF/expsigma/De",
                                                          "wTOF/expsigma/Tr", "wTOF/expsigma/He", "wTOF/expsigma/Al"};
  static constexpr std::string_view hnsigma_wTOF[Np] = {"wTOF/nsigma/El", "wTOF/nsigma/Mu", "wTOF/nsigma/Pi",
                                                        "wTOF/nsigma/Ka", "wTOF/nsigma/Pr", "wTOF/nsigma/De",
                                                        "wTOF/nsigma/Tr", "wTOF/nsigma/He", "wTOF/nsigma/Al"};
  static constexpr std::string_view hnsigma_pt_wTOF[Np] = {"wTOF/nsigma/pt/El", "wTOF/nsigma/pt/Mu", "wTOF/nsigma/pt/Pi",
                                                           "wTOF/nsigma/pt/Ka", "wTOF/nsigma/pt/Pr", "wTOF/nsigma/pt/De",
                                                           "wTOF/nsigma/pt/Tr", "wTOF/nsigma/pt/He", "wTOF/nsigma/pt/Al"};
  static constexpr std::string_view hnsigma_pt_pos_wTOF[Np] = {"wTOF/nsigma/pt/pos/El", "wTOF/nsigma/pt/pos/Mu", "wTOF/nsigma/pt/pos/Pi",
                                                               "wTOF/nsigma/pt/pos/Ka", "wTOF/nsigma/pt/pos/Pr", "wTOF/nsigma/pt/pos/De",
                                                               "wTOF/nsigma/pt/pos/Tr", "wTOF/nsigma/pt/pos/He", "wTOF/nsigma/pt/pos/Al"};
  static constexpr std::string_view hnsigma_pt_neg_wTOF[Np] = {"wTOF/nsigma/pt/neg/El", "wTOF/nsigma/pt/neg/Mu", "wTOF/nsigma/pt/neg/Pi",
                                                               "wTOF/nsigma/pt/neg/Ka", "wTOF/nsigma/pt/neg/Pr", "wTOF/nsigma/pt/neg/De",
                                                               "wTOF/nsigma/pt/neg/Tr", "wTOF/nsigma/pt/neg/He", "wTOF/nsigma/pt/neg/Al"};
  static constexpr std::string_view hsignal_wTOF[Np] = {"wTOF/signal/El", "wTOF/signal/Mu", "wTOF/signal/Pi",
                                                        "wTOF/signal/Ka", "wTOF/signal/Pr", "wTOF/signal/De",
                                                        "wTOF/signal/Tr", "wTOF/signal/He", "wTOF/signal/Al"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<int> nBinsDelta{"nBinsDelta", 200, "Number of bins for the Delta"};
  Configurable<float> minDelta{"minDelta", -1000.f, "Minimum Delta in range"};
  Configurable<float> maxDelta{"maxDelta", 1000.f, "Maximum Delta in range"};
  Configurable<int> nBinsExpSigma{"nBinsExpSigma", 200, "Number of bins for the ExpSigma"};
  Configurable<float> minExpSigma{"minExpSigma", 0.f, "Minimum ExpSigma in range"};
  Configurable<float> maxExpSigma{"maxExpSigma", 200.f, "Maximum ExpSigma in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Maximum NSigma in range"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<bool> applyTrackCut{"applyTrackCut", false, "Flag to apply standard track cuts"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};

  template <o2::track::PID::ID id>
  void initPerParticle(const AxisSpec& pAxis, const AxisSpec& ptAxis)
  {
    static_assert(id >= 0 && id <= PID::Alpha && "Particle index outside limits");
    bool enableTOFHistos = false;
    bool enableFullHistos = false;
    int enabledProcesses = 0;
    switch (id) { // Skipping disabled particles
#define particleCase(particleId)                                                                     \
  case PID::particleId:                                                                              \
    if (!doprocess##particleId && !doprocessFull##particleId && !doprocessFullWithTOF##particleId) { \
      return;                                                                                        \
    }                                                                                                \
    if (doprocess##particleId) {                                                                     \
      enabledProcesses++;                                                                            \
    }                                                                                                \
    if (doprocessFull##particleId) {                                                                 \
      enableFullHistos = true;                                                                       \
      enabledProcesses++;                                                                            \
    }                                                                                                \
    if (doprocessFullWithTOF##particleId) {                                                          \
      enableFullHistos = true;                                                                       \
      enableTOFHistos = true;                                                                        \
      enabledProcesses++;                                                                            \
    }                                                                                                \
    LOGF(info, "Enabled TPC QA for %s %s", #particleId, pT[id]);                                     \
    break;

      particleCase(Electron);
      particleCase(Muon);
      particleCase(Pion);
      particleCase(Kaon);
      particleCase(Proton);
      particleCase(Deuteron);
      particleCase(Triton);
      particleCase(Helium3);
      particleCase(Alpha);
#undef particleCase
    }
    if (enabledProcesses != 1) {
      LOG(fatal) << "Cannot enable more than one process function per particle, check and retry!";
    }

    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TPC}(%s)", pT[id]);
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
    histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigma_pt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigma_pt_pos[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigma_pt_neg[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});

    if (!enableFullHistos) { // Enabling only NSigma for tiny tables
      return;
    }

    // Exp signal
    const AxisSpec expAxis{1000, 0, 1000, Form("d#it{E}/d#it{x}_(%s) A.U.", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)", pT[id])};
    axisTitle = Form("#Delta^{TPC}(%s)", pT[id]);
    histos.add(hdelta[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
    histos.add(hdelta_pt_pos[id].data(), Form("%s Positive tracks", axisTitle), kTH2F, {ptAxis, deltaAxis});
    histos.add(hdelta_pt_neg[id].data(), Form("%s Negative tracks", axisTitle), kTH2F, {ptAxis, deltaAxis});

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TPC}(%s)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});

    if (!enableTOFHistos) { // Returning if the plots with TOF are not requested
      return;
    }

    histos.add(hexpected_wTOF[id].data(), "With TOF", kTH2F, {pAxis, expAxis});
    histos.add(hdelta_wTOF[id].data(), "With TOF", kTH2F, {pAxis, deltaAxis});
    histos.add(hdelta_pt_pos_wTOF[id].data(), "With TOF Positive", kTH2F, {ptAxis, deltaAxis});
    histos.add(hdelta_pt_neg_wTOF[id].data(), "With TOF Negative", kTH2F, {ptAxis, deltaAxis});
    histos.add(hexpsigma_wTOF[id].data(), "With TOF", kTH2F, {pAxis, expSigmaAxis});
    histos.add(hnsigma_wTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigma_pt_wTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigma_pt_pos_wTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigma_pt_neg_wTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {ptAxis, nSigmaAxis});

    const AxisSpec dedxAxis{1000, 0, 1000, "d#it{E}/d#it{x} A.U."};
    histos.add(hsignal_wTOF[id].data(), "With TOF", kTH2F, {pAxis, dedxAxis});
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec multAxis{1000, 0.f, 1000.f, "Track multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec phiAxis{100, 0, TMath::TwoPi(), "#it{#phi}"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec pAxisPosNeg{2 * nBinsP, -maxP, maxP, "#it{p}/z (GeV/#it{c})"};
    const AxisSpec ptAxisPosNeg{2 * nBinsP, -maxP, maxP, "#it{p}_{T}/z (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
    }
    const AxisSpec dedxAxis{5000, 0, 5000, "d#it{E}/d#it{x} A.U."};

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1F, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed mult.");
    h->GetXaxis()->SetBinLabel(4, "Passed vtx Z");

    h = histos.add<TH1>("event/trackselection", "", kTH1F, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(3, "hasITS");
    h->GetXaxis()->SetBinLabel(4, "hasTPC");

    histos.add("event/vertexz", "", kTH1F, {vtxZAxis});
    h = histos.add<TH1>("event/particlehypo", "", kTH1F, {{10, 0, 10, "PID in tracking"}});
    for (int i = 0; i < 9; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PID::getName(i));
    }
    histos.add("event/trackmultiplicity", "", kTH1F, {multAxis});
    histos.add("event/tpcsignal", "", kTH2F, {pAxis, dedxAxis});
    histos.add("event/signedtpcsignal", "", kTH2F, {pAxisPosNeg, dedxAxis});
    histos.add("event/tpcsignalvspt", "", kTH2F, {ptAxis, dedxAxis});
    histos.add("event/signedtpcsignalvspt", "", kTH2F, {ptAxisPosNeg, dedxAxis});
    histos.add("event/eta", "", kTH1F, {etaAxis});
    histos.add("event/phi", "", kTH1F, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1F, {lAxis});
    histos.add("event/pt", "", kTH1F, {ptAxis});
    histos.add("event/p", "", kTH1F, {pAxis});

    static_for<0, 8>([&](auto i) {
      initPerParticle<i>(pAxis, ptAxis);
    });
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& tracks)
  {

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 1);
    }
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 2);
    }

    // Computing Multiplicity first
    int ntracks = 0;
    if constexpr (fillHistograms) {
      for (auto t : tracks) {
        if (applyTrackCut && !t.isGlobalTrack()) {
          continue;
        }
        ntracks += 1;
      }
      histos.fill(HIST("event/evsel"), 3);
    }
    if (abs(collision.posZ()) > 10.f) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 4);
      histos.fill(HIST("event/vertexz"), collision.posZ());
      histos.fill(HIST("event/trackmultiplicity"), ntracks);
    }
    return true;
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isTrackSelected(const CollisionType& collision, const TrackType& track)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 1.f);
    }
    if (!track.isGlobalTrack()) { // Skipping non global tracks
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 2.f);
    }
    if (!track.hasITS()) { // Skipping tracks without ITS
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 3.f);
    }
    if (!track.hasTPC()) { // Skipping tracks without TPC
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 4.f);
      histos.fill(HIST("event/particlehypo"), track.pidForTracking());
      histos.fill(HIST("event/tpcsignal"), track.tpcInnerParam(), track.tpcSignal());
      histos.fill(HIST("event/signedtpcsignal"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
      histos.fill(HIST("event/tpcsignalvspt"), track.pt(), track.tpcSignal());
      histos.fill(HIST("event/signedtpcsignalvspt"), track.pt() * track.sign(), track.tpcSignal());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/phi"), track.phi());
      histos.fill(HIST("event/etaphi"), track.eta(), track.phi());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
      // histos.fill(HIST("event/ptreso"), track.p(), track.sigma1Pt() * track.pt() * track.pt());
    }
    return true;
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  void process(CollisionCandidate const& collision,
               soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::TrackSelection> const& tracks)
  {
    isEventSelected<true>(collision, tracks);
    for (auto t : tracks) {
      isTrackSelected<true>(collision, t);
    }
  }

  template <o2::track::PID::ID id, bool fillFullHistograms,
            bool fillWithTOFHistograms,
            typename TrackType>
  void processSingleParticle(CollisionCandidate const& collision,
                             TrackType const& tracks)
  {
    if (!isEventSelected<false>(collision, tracks)) {
      return;
    }

    for (auto t : tracks) {
      if (!isTrackSelected<false>(collision, t)) {
        continue;
      }

      if (applyRapidityCut) {
        if (abs(t.rapidity(PID::getMass(id))) > 0.5) {
          continue;
        }
      }

      const auto nsigma = o2::aod::pidutils::tpcNSigma<id>(t);
      histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
      histos.fill(HIST(hnsigma_pt[id]), t.pt(), nsigma);
      if (t.sign() > 0) {
        histos.fill(HIST(hnsigma_pt_pos[id]), t.pt(), nsigma);
      } else {
        histos.fill(HIST(hnsigma_pt_neg[id]), t.pt(), nsigma);
      }

      if constexpr (fillFullHistograms) {
        const auto& diff = o2::aod::pidutils::tpcExpSignalDiff<id>(t);
        // Fill histograms
        histos.fill(HIST(hexpected[id]), t.tpcInnerParam(), t.tpcSignal() - diff);
        histos.fill(HIST(hdelta[id]), t.tpcInnerParam(), diff);
        if (t.sign() > 0) {
          histos.fill(HIST(hdelta_pt_pos[id]), t.pt(), diff);
        } else {
          histos.fill(HIST(hdelta_pt_neg[id]), t.pt(), diff);
        }
        histos.fill(HIST(hexpsigma[id]), t.tpcInnerParam(), o2::aod::pidutils::tpcExpSigma<id>(t));
        if constexpr (fillWithTOFHistograms) {
          if (std::abs(o2::aod::pidutils::tofNSigma<id>(t)) < 3.f) {
            histos.fill(HIST(hexpected_wTOF[id]), t.tpcInnerParam(), t.tpcSignal() - diff);
            histos.fill(HIST(hdelta_wTOF[id]), t.tpcInnerParam(), diff);
            histos.fill(HIST(hexpsigma_wTOF[id]), t.p(), o2::aod::pidutils::tpcExpSigma<id>(t));
          }
        }
      }
      if constexpr (fillWithTOFHistograms) { // Filling nsigma (common to full and tiny)
        const auto& nsigmatof = o2::aod::pidutils::tofNSigma<id>(t);
        if (std::abs(nsigmatof) < 3.f) {
          histos.fill(HIST(hnsigma_wTOF[id]), t.p(), nsigma);
          histos.fill(HIST(hnsigma_pt_wTOF[id]), t.pt(), nsigma);
          histos.fill(HIST(hsignal_wTOF[id]), t.tpcInnerParam(), t.tpcSignal());
          // histos.fill(HIST("event/signedtpcsignal"), t.tpcInnerParam() * t.sign(), t.tpcSignal());
        }
      }
    }
  }

  // QA of nsigma only tables
#define makeProcessFunction(inputPid, particleId)                                        \
  void process##particleId(CollisionCandidate const& collision,                          \
                           soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, \
                                     inputPid> const& tracks)                            \
  {                                                                                      \
    processSingleParticle<PID::particleId, false, false>(collision, tracks);             \
  }                                                                                      \
  PROCESS_SWITCH(tpcPidQa, process##particleId, Form("Process for the %s hypothesis for TPC NSigma QA", #particleId), false);

  makeProcessFunction(aod::pidTPCEl, Electron);
  makeProcessFunction(aod::pidTPCMu, Muon);
  makeProcessFunction(aod::pidTPCPi, Pion);
  makeProcessFunction(aod::pidTPCKa, Kaon);
  makeProcessFunction(aod::pidTPCPr, Proton);
  makeProcessFunction(aod::pidTPCDe, Deuteron);
  makeProcessFunction(aod::pidTPCTr, Triton);
  makeProcessFunction(aod::pidTPCHe, Helium3);
  makeProcessFunction(aod::pidTPCAl, Alpha);
#undef makeProcessFunction

// QA of full tables
#define makeProcessFunction(inputPid, particleId)                                            \
  void processFull##particleId(CollisionCandidate const& collision,                          \
                               soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, \
                                         inputPid> const& tracks)                            \
  {                                                                                          \
    processSingleParticle<PID::particleId, true, false>(collision, tracks);                  \
  }                                                                                          \
  PROCESS_SWITCH(tpcPidQa, processFull##particleId, Form("Process for the %s hypothesis for full TPC PID QA", #particleId), false);

  makeProcessFunction(aod::pidTPCFullEl, Electron);
  makeProcessFunction(aod::pidTPCFullMu, Muon);
  makeProcessFunction(aod::pidTPCFullPi, Pion);
  makeProcessFunction(aod::pidTPCFullKa, Kaon);
  makeProcessFunction(aod::pidTPCFullPr, Proton);
  makeProcessFunction(aod::pidTPCFullDe, Deuteron);
  makeProcessFunction(aod::pidTPCFullTr, Triton);
  makeProcessFunction(aod::pidTPCFullHe, Helium3);
  makeProcessFunction(aod::pidTPCFullAl, Alpha);
#undef makeProcessFunction

  // QA of full tables with TOF information
#define makeProcessFunction(inputPid, inputPidTOF, particleId)                                                                            \
  void processFullWithTOF##particleId(CollisionCandidate const& collision,                                                                \
                                      soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, inputPid, inputPidTOF> const& tracks) \
  {                                                                                                                                       \
    processSingleParticle<PID::particleId, true, true>(collision, tracks);                                                                \
  }                                                                                                                                       \
  PROCESS_SWITCH(tpcPidQa, processFullWithTOF##particleId, Form("Process for the %s hypothesis for full TPC PID QA", #particleId), false);

  makeProcessFunction(aod::pidTPCFullEl, aod::pidTOFFullEl, Electron);
  makeProcessFunction(aod::pidTPCFullMu, aod::pidTOFFullMu, Muon);
  makeProcessFunction(aod::pidTPCFullPi, aod::pidTOFFullPi, Pion);
  makeProcessFunction(aod::pidTPCFullKa, aod::pidTOFFullKa, Kaon);
  makeProcessFunction(aod::pidTPCFullPr, aod::pidTOFFullPr, Proton);
  makeProcessFunction(aod::pidTPCFullDe, aod::pidTOFFullDe, Deuteron);
  makeProcessFunction(aod::pidTPCFullTr, aod::pidTOFFullTr, Triton);
  makeProcessFunction(aod::pidTPCFullHe, aod::pidTOFFullHe, Helium3);
  makeProcessFunction(aod::pidTPCFullAl, aod::pidTOFFullAl, Alpha);
#undef makeProcessFunction
};
