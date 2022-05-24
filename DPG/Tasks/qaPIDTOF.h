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
  Configurable<bool> enableEvTimeSplitting{"enableEvTimeSplitting", false, "Flag to enable histograms splitting depending on the Event Time used"};

  template <o2::track::PID::ID id>
  void initPerParticle(const AxisSpec& pAxis, const AxisSpec& ptAxis)
  {
    static_assert(id >= 0 && id <= PID::Alpha && "Particle index outside limits");
    switch (id) { // Skipping disabled particles
#define particleCase(particleId)                                 \
  case PID::particleId:                                          \
    if (!doprocess##particleId && !doprocessFull##particleId) {  \
      return;                                                    \
    }                                                            \
    LOGF(info, "Enabled TOF QA for %s %s", #particleId, pT[id]); \
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

    auto addHistogram = [&](const auto& name, const auto& title, const auto& xAxis, const auto& yAxis) {
      if (!enableEvTimeSplitting) {
        histos.add(name.data(), title, kTH2F, {xAxis, yAxis});
        return;
      }
      const AxisSpec nEvTimeTypeAxis{4, 0.5, 4.5, "Event Time used"};
      auto histo = histos.add<TH3>(name.data(), title, kTH3F, {xAxis, yAxis, nEvTimeTypeAxis});
      histo->GetZaxis()->SetBinLabel(1, "No Ev. Time");
      histo->GetZaxis()->SetBinLabel(2, "TOF");
      histo->GetZaxis()->SetBinLabel(3, "FT0");
      histo->GetZaxis()->SetBinLabel(4, "FT0+TOF");
    };

    // Exp signal
    const AxisSpec expAxis{1000, 0, 2e6, Form("t_{exp}(%s) (ps)", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("t-t_{ev}-t_{exp}(%s) (ps)", pT[id])};
    addHistogram(hexpected_diff[id], Form("#Delta^{TOF}(%s)", pT[id]), pAxis, deltaAxis);

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TOF}(%s) (ps)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});

    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TOF}(%s)", pT[id]);
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
    addHistogram(hnsigma[id], axisTitle, pAxis, nSigmaAxis);
    addHistogram(hnsigmapt[id], axisTitle, ptAxis, nSigmaAxis);
    addHistogram(hnsigmapospt[id], axisTitle, ptAxis, nSigmaAxis);
    addHistogram(hnsigmanegpt[id], axisTitle, ptAxis, nSigmaAxis);
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
    int evTimeIndex = 1; // Index for the event time type

    for (auto t : tracks) {
      if (!isTrackSelected<false>(collision, t)) {
        continue;
      }

      if (applyRapidityCut) {
        if (abs(t.rapidity(PID::getMass(id))) > 0.5) {
          continue;
        }
      }

      const auto nsigma = o2::aod::pidutils::tofNSigma<id>(t);
      if (!enableEvTimeSplitting) {
        histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
        histos.fill(HIST(hnsigmapt[id]), t.pt(), nsigma);
        if (t.sign() > 0) {
          histos.fill(HIST(hnsigmapospt[id]), t.pt(), nsigma);
        } else {
          histos.fill(HIST(hnsigmanegpt[id]), t.pt(), nsigma);
        }
      } else {
        evTimeIndex = 1;
        if (t.isEvTimeTOF() && t.isEvTimeT0AC()) {
          evTimeIndex = 4;
        } else if (t.isEvTimeT0AC()) {
          evTimeIndex = 3;
        } else if (t.isEvTimeTOF()) {
          evTimeIndex = 2;
        }

        histos.fill(HIST(hnsigma[id]), t.p(), nsigma, evTimeIndex);
        histos.fill(HIST(hnsigmapt[id]), t.pt(), nsigma, evTimeIndex);
        if (t.sign() > 0) {
          histos.fill(HIST(hnsigmapospt[id]), t.pt(), nsigma, evTimeIndex);
        } else {
          histos.fill(HIST(hnsigmanegpt[id]), t.pt(), nsigma, evTimeIndex);
        }
      }
      if constexpr (fillFullHistograms) {
        const float tof = t.tofSignal() - collisionTime_ps;
        const auto diff = o2::aod::pidutils::tofExpSignalDiff<id>(t);
        histos.fill(HIST(hexpected[id]), t.p(), tof - diff);
        if (!enableEvTimeSplitting) {
          histos.fill(HIST(hexpected_diff[id]), t.p(), diff);
        } else {
          histos.fill(HIST(hexpected_diff[id]), t.p(), diff, evTimeIndex);
        }
        histos.fill(HIST(hexpsigma[id]), t.p(), o2::aod::pidutils::tofExpSigma<id>(t));
      }
    }
  }

  // QA of nsigma only tables
#define makeProcessFunction(inputPid, particleId)                                                  \
  void process##particleId(CollisionCandidate const& collision,                                    \
                           soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,           \
                                     aod::pidEvTimeFlags, aod::TOFSignal, inputPid> const& tracks) \
  {                                                                                                \
    processSingleParticle<PID::particleId, false>(collision, tracks);                              \
  }                                                                                                \
  PROCESS_SWITCH(tofPidQa, process##particleId, Form("Process for the %s hypothesis for TOF NSigma QA", #particleId), false);

  makeProcessFunction(aod::pidTOFEl, Electron);
  makeProcessFunction(aod::pidTOFMu, Muon);
  makeProcessFunction(aod::pidTOFPi, Pion);
  makeProcessFunction(aod::pidTOFKa, Kaon);
  makeProcessFunction(aod::pidTOFPr, Proton);
  makeProcessFunction(aod::pidTOFDe, Deuteron);
  makeProcessFunction(aod::pidTOFTr, Triton);
  makeProcessFunction(aod::pidTOFHe, Helium3);
  makeProcessFunction(aod::pidTOFAl, Alpha);
#undef makeProcessFunction

// QA of full tables
#define makeProcessFunction(inputPid, particleId)                                                      \
  void processFull##particleId(CollisionCandidate const& collision,                                    \
                               soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,           \
                                         aod::pidEvTimeFlags, aod::TOFSignal, inputPid> const& tracks) \
  {                                                                                                    \
    processSingleParticle<PID::particleId, true>(collision, tracks);                                   \
  }                                                                                                    \
  PROCESS_SWITCH(tofPidQa, processFull##particleId, Form("Process for the %s hypothesis for full TOF PID QA", #particleId), false);

  makeProcessFunction(aod::pidTOFFullEl, Electron);
  makeProcessFunction(aod::pidTOFFullMu, Muon);
  makeProcessFunction(aod::pidTOFFullPi, Pion);
  makeProcessFunction(aod::pidTOFFullKa, Kaon);
  makeProcessFunction(aod::pidTOFFullPr, Proton);
  makeProcessFunction(aod::pidTOFFullDe, Deuteron);
  makeProcessFunction(aod::pidTOFFullTr, Triton);
  makeProcessFunction(aod::pidTOFFullHe, Helium3);
  makeProcessFunction(aod::pidTOFFullAl, Alpha);
#undef makeProcessFunction
};

/// Task to produce the TOF QA plots for Beta
struct tofPidBetaQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hexpected_diff[Np] = {"expected_diff/El", "expected_diff/Mu", "expected_diff/Pi",
                                                          "expected_diff/Ka", "expected_diff/Pr", "expected_diff/De",
                                                          "expected_diff/Tr", "expected_diff/He", "expected_diff/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};
  Configurable<int> nBinsBeta{"nBinsBeta", 4000, "Number of bins for the beta"};
  Configurable<float> minBeta{"minBeta", 0, "Minimum beta in range"};
  Configurable<float> maxBeta{"maxBeta", 2.f, "Maximum beta in range"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal"};
    const AxisSpec betaAxis{nBinsBeta, minBeta, maxBeta, "TOF #beta"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogaritmic();
      pAxis.makeLogaritmic();
    }

    // Event properties
    histos.add("event/tofsignal", "", HistType::kTH2F, {pAxis, tofAxis});
    histos.add("event/tofbeta", "", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeTOF", "", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/eta", "", HistType::kTH1F, {etaAxis});
    histos.add("event/length", "", HistType::kTH1F, {lAxis});
    histos.add("event/pt", "", HistType::kTH1F, {ptAxis});
    histos.add("event/p", "", HistType::kTH1F, {pAxis});
  }

  template <uint8_t i, typename T>
  void fillParticleHistos(const T& t, const float tof, const float exp_diff, const float nsigma)
  {
    histos.fill(HIST(hexpected[i]), t.p(), tof - exp_diff);
    histos.fill(HIST(hexpected_diff[i]), t.p(), exp_diff);
    histos.fill(HIST(hnsigma[i]), t.p(), nsigma);
  }
  void process(soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::pidEvTimeFlags> const& tracks,
               aod::Collisions const&)
  {
    for (auto const& track : tracks) {

      if (!track.hasTOF()) { // Skipping tracks without TOF
        continue;
      }
      if (!track.isGlobalTrack()) {
        continue;
      }
      if (track.isEvTimeTOF()) {
        histos.fill(HIST("event/tofbetaEvTimeTOF"), track.p(), track.beta());
      }
      histos.fill(HIST("event/tofbeta"), track.p(), track.beta());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
    }
  }
};
