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
          return;
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

  template <typename pidhypothesis>
  using TrackCandidate = soa::Join<aod::Tracks, aod::TracksExtra, pidhypothesis, aod::TrackSelection>;

  // QA of nsigma only tables
  void processElectron(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTPCEl> const& tracks)
  {
    processSingleParticle<PID::Electron, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processElectron, "Process for the Electron hypothesis", false);

  void processMuon(CollisionCandidate const& collision,
                   TrackCandidate<aod::pidTPCMu> const& tracks)
  {
    processSingleParticle<PID::Muon, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processMuon, "Process for the Muon hypothesis", false);

  void processPion(CollisionCandidate const& collision,
                   TrackCandidate<aod::pidTPCPi> const& tracks)
  {
    processSingleParticle<PID::Pion, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processPion, "Process for the Pion hypothesis", false);

  void processKaon(CollisionCandidate const& collision,
                   TrackCandidate<aod::pidTPCKa> const& tracks)
  {
    processSingleParticle<PID::Kaon, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processKaon, "Process for the Kaon hypothesis", false);

  void processProton(CollisionCandidate const& collision,
                     TrackCandidate<aod::pidTPCPr> const& tracks)
  {
    processSingleParticle<PID::Proton, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processProton, "Process for the Proton hypothesis", false);

  void processDeuteron(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTPCDe> const& tracks)
  {
    processSingleParticle<PID::Deuteron, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processDeuteron, "Process for the Deuteron hypothesis", false);

  void processTriton(CollisionCandidate const& collision,
                     TrackCandidate<aod::pidTPCTr> const& tracks)
  {
    processSingleParticle<PID::Triton, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processTriton, "Process for the Triton hypothesis", false);

  void processHelium3(CollisionCandidate const& collision,
                      TrackCandidate<aod::pidTPCHe> const& tracks)
  {
    processSingleParticle<PID::Helium3, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processHelium3, "Process for the Helium3 hypothesis", false);

  void processAlpha(CollisionCandidate const& collision,
                    TrackCandidate<aod::pidTPCAl> const& tracks)
  {
    processSingleParticle<PID::Alpha, false, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processAlpha, "Process for the Alpha hypothesis", false);

  // QA of full tables
  void processFullElectron(CollisionCandidate const& collision,
                           TrackCandidate<aod::pidTPCFullEl> const& tracks)
  {
    processSingleParticle<PID::Electron, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullElectron, "Process for the Electron hypothesis for full PID information", false);

  void processFullMuon(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTPCFullMu> const& tracks)
  {
    processSingleParticle<PID::Muon, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullMuon, "Process for the Muon hypothesis for full PID information", false);

  void processFullPion(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTPCFullPi> const& tracks)
  {
    processSingleParticle<PID::Pion, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullPion, "Process for the Pion hypothesis for full PID information", false);

  void processFullKaon(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTPCFullKa> const& tracks)
  {
    processSingleParticle<PID::Kaon, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullKaon, "Process for the Kaon hypothesis for full PID information", false);

  void processFullProton(CollisionCandidate const& collision,
                         TrackCandidate<aod::pidTPCFullPr> const& tracks)
  {
    processSingleParticle<PID::Proton, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullProton, "Process for the Proton hypothesis for full PID information", false);

  void processFullDeuteron(CollisionCandidate const& collision,
                           TrackCandidate<aod::pidTPCFullDe> const& tracks)
  {
    processSingleParticle<PID::Deuteron, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullDeuteron, "Process for the Deuteron hypothesis for full PID information", false);

  void processFullTriton(CollisionCandidate const& collision,
                         TrackCandidate<aod::pidTPCFullTr> const& tracks)
  {
    processSingleParticle<PID::Triton, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullTriton, "Process for the Triton hypothesis for full PID information", false);

  void processFullHelium3(CollisionCandidate const& collision,
                          TrackCandidate<aod::pidTPCFullHe> const& tracks)
  {
    processSingleParticle<PID::Helium3, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullHelium3, "Process for the Helium3 hypothesis for full PID information", false);

  void processFullAlpha(CollisionCandidate const& collision,
                        TrackCandidate<aod::pidTPCFullAl> const& tracks)
  {
    processSingleParticle<PID::Alpha, true, false>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullAlpha, "Process for the Alpha hypothesis for full PID information", false);

  // QA of full tables with TOF information
  void processFullWithTOFElectron(CollisionCandidate const& collision,
                                  TrackCandidate<soa::Join<aod::pidTOFFullEl, aod::pidTPCFullEl>> const& tracks)
  {
    processSingleParticle<PID::Electron, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFElectron, "Process for the Electron hypothesis for full PID information", false);

  void processFullWithTOFMuon(CollisionCandidate const& collision,
                              TrackCandidate<soa::Join<aod::pidTOFFullMu, aod::pidTPCFullMu>> const& tracks)
  {
    processSingleParticle<PID::Muon, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFMuon, "Process for the Muon hypothesis for full PID information", false);

  void processFullWithTOFPion(CollisionCandidate const& collision,
                              TrackCandidate<soa::Join<aod::pidTOFFullPi, aod::pidTPCFullPi>> const& tracks)
  {
    processSingleParticle<PID::Pion, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFPion, "Process for the Pion hypothesis for full PID information", false);

  void processFullWithTOFKaon(CollisionCandidate const& collision,
                              TrackCandidate<soa::Join<aod::pidTOFFullKa, aod::pidTPCFullKa>> const& tracks)
  {
    processSingleParticle<PID::Kaon, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFKaon, "Process for the Kaon hypothesis for full PID information", false);

  void processFullWithTOFProton(CollisionCandidate const& collision,
                                TrackCandidate<soa::Join<aod::pidTOFFullPr, aod::pidTPCFullPr>> const& tracks)
  {
    processSingleParticle<PID::Proton, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFProton, "Process for the Proton hypothesis for full PID information", false);

  void processFullWithTOFDeuteron(CollisionCandidate const& collision,
                                  TrackCandidate<soa::Join<aod::pidTOFFullDe, aod::pidTPCFullDe>> const& tracks)
  {
    processSingleParticle<PID::Deuteron, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFDeuteron, "Process for the Deuteron hypothesis for full PID information", false);

  void processFullWithTOFTriton(CollisionCandidate const& collision,
                                TrackCandidate<soa::Join<aod::pidTOFFullTr, aod::pidTPCFullTr>> const& tracks)
  {
    processSingleParticle<PID::Triton, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFTriton, "Process for the Triton hypothesis for full PID information", false);

  void processFullWithTOFHelium3(CollisionCandidate const& collision,
                                 TrackCandidate<soa::Join<aod::pidTOFFullHe, aod::pidTPCFullHe>> const& tracks)
  {
    processSingleParticle<PID::Helium3, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFHelium3, "Process for the Helium3 hypothesis for full PID information", false);

  void processFullWithTOFAlpha(CollisionCandidate const& collision,
                               TrackCandidate<soa::Join<aod::pidTOFFullAl, aod::pidTPCFullAl>> const& tracks)
  {
    processSingleParticle<PID::Alpha, true, true>(collision, tracks);
  }
  PROCESS_SWITCH(tpcPidQa, processFullWithTOFAlpha, "Process for the Alpha hypothesis for full PID information", false);
};
