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
/// \file   qaTOF.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Header file for QA tasks of the TOF PID quantities
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

/// Task to produce the TOF QA plots
struct tofPidQa {
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

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};
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
  void initPerParticle(const AxisSpec& pAxis, const AxisSpec& ptAxis)
  {
    // Exp signal
    const AxisSpec expAxis{1000, 0, 2e6, Form("t_{exp}(%s) (ps)", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("t-t_{ev}-t_{exp}(%s) (ps)", pT[id])};
    histos.add(hexpected_diff[id].data(), "", kTH2F, {pAxis, deltaAxis});

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TOF}(%s) (ps)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});

    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TOF}(%s)", pT[id]);
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
    histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigmapt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmapospt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmanegpt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec multAxis{100, 0, 100, "TOF multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal (ps)"};
    const AxisSpec etaAxis{100, -1, 1, "#it{#eta}"};
    const AxisSpec phiAxis{100, 0, TMath::TwoPi(), "#it{#phi}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec colTimeResoAxis{100, 0, 1000, "#sigma_{Collision time} (ps)"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec pExpAxis{nBinsP, minP, maxP, "#it{p}_{Exp. TOF} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogaritmic();
      pAxis.makeLogaritmic();
      pExpAxis.makeLogaritmic();
    }

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
    h->GetXaxis()->SetBinLabel(5, "hasTOF");

    histos.add("event/vertexz", "", kTH1F, {vtxZAxis});
    h = histos.add<TH1>("event/particlehypo", "", kTH1F, {{10, 0, 10, "PID in tracking"}});
    for (int i = 0; i < 9; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PID::getName(i));
    }
    histos.add("event/trackmultiplicity", "", kTH1F, {multAxis});
    histos.add("event/tofmultiplicity", "", kTH1F, {multAxis});
    histos.add("event/colltime", "", kTH1F, {colTimeAxis});
    histos.add("event/colltimereso", "", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/tofsignal", "", kTH2F, {pAxis, tofAxis});
    histos.add("event/pexp", "", kTH2F, {pAxis, pExpAxis});
    histos.add("event/eta", "", kTH1F, {etaAxis});
    histos.add("event/phi", "", kTH1F, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1F, {lAxis});
    histos.add("event/pt", "", kTH1F, {ptAxis});
    histos.add("event/p", "", kTH1F, {pAxis});
    // histos.add("event/ptreso", "", kTH2F, {pAxis, ptResoAxis});

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
    float ntracks = 0;
    int tofmult = 0;
    if constexpr (fillHistograms) {
      for (auto t : tracks) {
        if (applyTrackCut && !t.isGlobalTrack()) {
          continue;
        }
        ntracks += 1;
        if (!t.hasTOF()) { // Skipping tracks without TOF
          continue;
        }
        tofmult++;
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
      histos.fill(HIST("event/tofmultiplicity"), tofmult);

      const float collisionTime_ps = collision.collisionTime() * 1000.f;
      histos.fill(HIST("event/colltime"), collisionTime_ps);
      histos.fill(HIST("event/colltimereso"), tofmult, collision.collisionTimeRes() * 1000.f);
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
    }
    if (!track.hasTOF()) { // Skipping tracks without TOF
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 5.f);
      histos.fill(HIST("event/particlehypo"), track.pidForTracking());
      histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      histos.fill(HIST("event/pexp"), track.p(), track.tofExpMom());
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
               soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection> const& tracks)
  {
    isEventSelected<true>(collision, tracks);
    for (auto t : tracks) {
      isTrackSelected<true>(collision, t);
    }
  }

  template <o2::track::PID::ID id, bool fillFullHistograms, typename TrackType>
  void processSingleParticle(CollisionCandidate const& collision,
                             TrackType const& tracks)
  {
    if (!isEventSelected<false>(collision, tracks)) {
      return;
    }
    const float collisionTime_ps = collision.collisionTime() * 1000.f;

    for (auto t : tracks) {
      if (!isTrackSelected<false>(collision, t)) {
        continue;
      }

      if (applyRapidityCut) {
        if (abs(t.rapidity(PID::getMass(id))) > 0.5) {
          return;
        }
      }

      const auto nsigma = o2::aod::pidutils::tofNSigma<id>(t);
      histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
      histos.fill(HIST(hnsigmapt[id]), t.pt(), nsigma);
      if (t.sign() > 0) {
        histos.fill(HIST(hnsigmapospt[id]), t.pt(), nsigma);
      } else {
        histos.fill(HIST(hnsigmanegpt[id]), t.pt(), nsigma);
      }
      if constexpr (fillFullHistograms) {
        const float tof = t.tofSignal() - collisionTime_ps;
        const auto diff = o2::aod::pidutils::tofExpSignalDiff<id>(t);
        histos.fill(HIST(hexpected[id]), t.p(), tof - diff);
        histos.fill(HIST(hexpected_diff[id]), t.p(), diff);
        histos.fill(HIST(hexpsigma[id]), t.p(), o2::aod::pidutils::tofExpSigma<id>(t));
      }
    }
  }

  template <typename pidhypothesis>
  using TrackCandidate = soa::Join<aod::Tracks, aod::TracksExtra, pidhypothesis, aod::TOFSignal, aod::TrackSelection>;

  // QA of nsigma only tables
  void processElectron(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTOFEl> const& tracks)
  {
    processSingleParticle<PID::Electron, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processElectron, "Process for the Electron hypothesis", false);

  void processMuon(CollisionCandidate const& collision,
                   TrackCandidate<aod::pidTOFMu> const& tracks)
  {
    processSingleParticle<PID::Muon, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processMuon, "Process for the Muon hypothesis", false);

  void processPion(CollisionCandidate const& collision,
                   TrackCandidate<aod::pidTOFPi> const& tracks)
  {
    processSingleParticle<PID::Pion, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processPion, "Process for the Pion hypothesis", false);

  void processKaon(CollisionCandidate const& collision,
                   TrackCandidate<aod::pidTOFKa> const& tracks)
  {
    processSingleParticle<PID::Kaon, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processKaon, "Process for the Kaon hypothesis", false);

  void processProton(CollisionCandidate const& collision,
                     TrackCandidate<aod::pidTOFPr> const& tracks)
  {
    processSingleParticle<PID::Proton, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processProton, "Process for the Proton hypothesis", false);

  void processDeuteron(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTOFDe> const& tracks)
  {
    processSingleParticle<PID::Deuteron, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processDeuteron, "Process for the Deuteron hypothesis", false);

  void processTriton(CollisionCandidate const& collision,
                     TrackCandidate<aod::pidTOFTr> const& tracks)
  {
    processSingleParticle<PID::Triton, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processTriton, "Process for the Triton hypothesis", false);

  void processHelium3(CollisionCandidate const& collision,
                      TrackCandidate<aod::pidTOFHe> const& tracks)
  {
    processSingleParticle<PID::Helium3, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processHelium3, "Process for the Helium3 hypothesis", false);

  void processAlpha(CollisionCandidate const& collision,
                    TrackCandidate<aod::pidTOFAl> const& tracks)
  {
    processSingleParticle<PID::Alpha, false>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processAlpha, "Process for the Alpha hypothesis", false);

  // QA of full tables
  void processFullElectron(CollisionCandidate const& collision,
                           TrackCandidate<aod::pidTOFFullEl> const& tracks)
  {
    processSingleParticle<PID::Electron, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullElectron, "Process for the Electron hypothesis", false);

  void processFullMuon(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTOFFullMu> const& tracks)
  {
    processSingleParticle<PID::Muon, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullMuon, "Process for the Muon hypothesis", false);

  void processFullPion(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTOFFullPi> const& tracks)
  {
    processSingleParticle<PID::Pion, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullPion, "Process for the Pion hypothesis", false);

  void processFullKaon(CollisionCandidate const& collision,
                       TrackCandidate<aod::pidTOFFullKa> const& tracks)
  {
    processSingleParticle<PID::Kaon, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullKaon, "Process for the Kaon hypothesis", false);

  void processFullProton(CollisionCandidate const& collision,
                         TrackCandidate<aod::pidTOFFullPr> const& tracks)
  {
    processSingleParticle<PID::Proton, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullProton, "Process for the Proton hypothesis", false);

  void processFullDeuteron(CollisionCandidate const& collision,
                           TrackCandidate<aod::pidTOFFullDe> const& tracks)
  {
    processSingleParticle<PID::Deuteron, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullDeuteron, "Process for the Deuteron hypothesis", false);

  void processFullTriton(CollisionCandidate const& collision,
                         TrackCandidate<aod::pidTOFFullTr> const& tracks)
  {
    processSingleParticle<PID::Triton, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullTriton, "Process for the Triton hypothesis", false);

  void processFullHelium3(CollisionCandidate const& collision,
                          TrackCandidate<aod::pidTOFFullHe> const& tracks)
  {
    processSingleParticle<PID::Helium3, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullHelium3, "Process for the Helium3 hypothesis", false);

  void processFullAlpha(CollisionCandidate const& collision,
                        TrackCandidate<aod::pidTOFFullAl> const& tracks)
  {
    processSingleParticle<PID::Alpha, true>(collision, tracks);
  }
  PROCESS_SWITCH(tofPidQa, processFullAlpha, "Process for the Alpha hypothesis", false);
};