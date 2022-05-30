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
/// \file   qaTPC.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Header file for QA tasks of the TPC PID quantities
///

#include "Common/Core/PID/PIDResponse.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

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
  static constexpr std::string_view hexpected_diff[Np] = {"expected_diff/El", "expected_diff/Mu", "expected_diff/Pi",
                                                          "expected_diff/Ka", "expected_diff/Pr", "expected_diff/De",
                                                          "expected_diff/Tr", "expected_diff/He", "expected_diff/Al"};
  static constexpr std::string_view hexpsigma[Np] = {"expsigma/El", "expsigma/Mu", "expsigma/Pi",
                                                     "expsigma/Ka", "expsigma/Pr", "expsigma/De",
                                                     "expsigma/Tr", "expsigma/He", "expsigma/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigmapt[Np] = {"nsigmapt/El", "nsigmapt/Mu", "nsigmapt/Pi",
                                                     "nsigmapt/Ka", "nsigmapt/Pr", "nsigmapt/De",
                                                     "nsigmapt/Tr", "nsigmapt/He", "nsigmapt/Al"};
  static constexpr std::string_view hnsigmapospt[Np] = {"nsigmapospt/El", "nsigmapospt/Mu", "nsigmapospt/Pi",
                                                        "nsigmapospt/Ka", "nsigmapospt/Pr", "nsigmapospt/De",
                                                        "nsigmapospt/Tr", "nsigmapospt/He", "nsigmapospt/Al"};
  static constexpr std::string_view hnsigmanegpt[Np] = {"nsigmanegpt/El", "nsigmanegpt/Mu", "nsigmanegpt/Pi",
                                                        "nsigmanegpt/Ka", "nsigmanegpt/Pr", "nsigmanegpt/De",
                                                        "nsigmanegpt/Tr", "nsigmanegpt/He", "nsigmanegpt/Al"};
  static constexpr std::string_view hsignal[Np] = {"signal/El", "signal/Mu", "signal/Pi",
                                                   "signal/Ka", "signal/Pr", "signal/De",
                                                   "signal/Tr", "signal/He", "signal/Al"};

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
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<bool> applyTrackCut{"applyTrackCut", false, "Flag to apply standard track cuts"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};

  template <o2::track::PID::ID id>
  void addParticleHistos(const AxisSpec& pAxis, const AxisSpec& ptAxis)
  {
    static_assert(id >= 0 && id <= PID::Alpha && "Particle index outside limits");
    switch (id) { // Skipping disabled particles
#define particleCase(particleId)                                                                     \
  case PID::particleId:                                                                              \
    if (!doprocess##particleId && !doprocessFull##particleId && !doprocessFullWithTOF##particleId) { \
      return;                                                                                        \
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
    // Exp signal
    const AxisSpec expAxis{1000, 0, 1000, Form("d#it{E}/d#it{x}_(%s) A.U.", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)", pT[id])};
    histos.add(hexpected_diff[id].data(), "", kTH2F, {pAxis, deltaAxis});

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TPC}(%s)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});

    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TPC}(%s)", pT[id]);
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
    histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigmapt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmapospt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmanegpt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});

    // Signal with TOF preselection
    const AxisSpec dedxAxis{1000, 0, 1000, "d#it{E}/d#it{x} A.U."};
    histos.add(hsignal[id].data(), "", kTH2F, {pAxis, dedxAxis});
  }

  void init(o2::framework::InitContext&)
  {

    const AxisSpec multAxis{1000, 0.f, 1000.f, "Track multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec phiAxis{100, 0, TMath::TwoPi(), "#it{#phi}"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec pAxisPosNeg{nBinsP, -maxP, maxP, "Signed #it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogaritmic();
      pAxis.makeLogaritmic();
    }
    const AxisSpec dedxAxis{1000, 0, 1000, "d#it{E}/d#it{x} A.U."};

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
    histos.add("event/eta", "", kTH1F, {etaAxis});
    histos.add("event/phi", "", kTH1F, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1F, {lAxis});
    histos.add("event/pt", "", kTH1F, {ptAxis});
    histos.add("event/p", "", kTH1F, {pAxis});

    static_for<0, 8>([&](auto i) {
      addParticleHistos<i>(pAxis, ptAxis);
    });
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& tracks)
  {

    histos.fill(HIST("event/evsel"), 1);
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }

    histos.fill(HIST("event/evsel"), 2);

    // Computing Multiplicity first
    float ntracks = 0;
    for (auto t : tracks) {
      if (applyTrackCut && !t.isGlobalTrack()) {
        continue;
      }
      ntracks += 1;
    }
    histos.fill(HIST("event/evsel"), 3);
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
               soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection> const& tracks)
  {
    isEventSelected<true>(collision, tracks);
    for (auto t : tracks) {
      isTrackSelected<true>(collision, t);
    }
  }

  template <o2::track::PID::ID id, bool fillFullHistograms, bool fillWithTOFHistograms, typename TrackType>
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

      const auto& nsigma = o2::aod::pidutils::tpcNSigma<id>(t);
      histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
      histos.fill(HIST(hnsigmapt[id]), t.pt(), nsigma);
      if (t.sign() > 0) {
        histos.fill(HIST(hnsigmapospt[id]), t.pt(), nsigma);
      } else {
        histos.fill(HIST(hnsigmanegpt[id]), t.pt(), nsigma);
      }

      if constexpr (fillFullHistograms) {
        const auto& diff = o2::aod::pidutils::tpcExpSignalDiff<id>(t);
        // Fill histograms
        histos.fill(HIST(hexpected[id]), t.tpcInnerParam(), t.tpcSignal() - diff);
        histos.fill(HIST(hexpected_diff[id]), t.tpcInnerParam(), diff);
        histos.fill(HIST(hexpsigma[id]), t.tpcInnerParam(), o2::aod::pidutils::tpcExpSigma<id>(t));
        if constexpr (fillWithTOFHistograms) {
          if (std::abs(o2::aod::pidutils::tofNSigma<id>(t)) < 3.f) {
            histos.fill(HIST(hexpected[id]), t.tpcInnerParam(), t.tpcSignal() - diff);
            histos.fill(HIST(hexpected_diff[id]), t.tpcInnerParam(), diff);
            histos.fill(HIST(hexpsigma[id]), t.p(), o2::aod::pidutils::tpcExpSigma<id>(t));
          }
        }
      }
      if constexpr (fillWithTOFHistograms) {
        const auto& nsigmatof = o2::aod::pidutils::tofNSigma<id>(t);
        if (std::abs(nsigmatof) < 3.f) {
          histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
          histos.fill(HIST(hnsigmapt[id]), t.pt(), nsigma);
          histos.fill(HIST(hsignal[id]), t.tpcInnerParam(), t.tpcSignal());
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
