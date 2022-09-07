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
  static constexpr std::string_view hexpected_diff[Np] = {"expected_diff/El", "expected_diff/Mu", "expected_diff/Pi",
                                                          "expected_diff/Ka", "expected_diff/Pr", "expected_diff/De",
                                                          "expected_diff/Tr", "expected_diff/He", "expected_diff/Al"};
  static constexpr std::string_view hexpected_diffptpos[Np] = {"expected_diffptpos/El", "expected_diffptpos/Mu", "expected_diffptpos/Pi",
                                                               "expected_diffptpos/Ka", "expected_diffptpos/Pr", "expected_diffptpos/De",
                                                               "expected_diffptpos/Tr", "expected_diffptpos/He", "expected_diffptpos/Al"};
  static constexpr std::string_view hexpected_diffptneg[Np] = {"expected_diffptneg/El", "expected_diffptneg/Mu", "expected_diffptneg/Pi",
                                                               "expected_diffptneg/Ka", "expected_diffptneg/Pr", "expected_diffptneg/De",
                                                               "expected_diffptneg/Tr", "expected_diffptneg/He", "expected_diffptneg/Al"};
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

  // With TOF
  static constexpr std::string_view hexpectedWithTOF[Np] = {"wTOF/expected/El", "wTOF/expected/Mu", "wTOF/expected/Pi",
                                                            "wTOF/expected/Ka", "wTOF/expected/Pr", "wTOF/expected/De",
                                                            "wTOF/expected/Tr", "wTOF/expected/He", "wTOF/expected/Al"};
  static constexpr std::string_view hexpected_diffWithTOF[Np] = {"wTOF/expected_diff/El", "wTOF/expected_diff/Mu", "wTOF/expected_diff/Pi",
                                                                 "wTOF/expected_diff/Ka", "wTOF/expected_diff/Pr", "wTOF/expected_diff/De",
                                                                 "wTOF/expected_diff/Tr", "wTOF/expected_diff/He", "wTOF/expected_diff/Al"};
  static constexpr std::string_view hexpected_diffptposWithTOF[Np] = {"wTOF/expected_diffptpos/El", "wTOF/expected_diffptpos/Mu", "wTOF/expected_diffptpos/Pi",
                                                                      "wTOF/expected_diffptpos/Ka", "wTOF/expected_diffptpos/Pr", "wTOF/expected_diffptpos/De",
                                                                      "wTOF/expected_diffptpos/Tr", "wTOF/expected_diffptpos/He", "wTOF/expected_diffptpos/Al"};
  static constexpr std::string_view hexpected_diffptnegWithTOF[Np] = {"wTOF/expected_diffptneg/El", "wTOF/expected_diffptneg/Mu", "wTOF/expected_diffptneg/Pi",
                                                                      "wTOF/expected_diffptneg/Ka", "wTOF/expected_diffptneg/Pr", "wTOF/expected_diffptneg/De",
                                                                      "wTOF/expected_diffptneg/Tr", "wTOF/expected_diffptneg/He", "wTOF/expected_diffptneg/Al"};
  static constexpr std::string_view hexpsigmaWithTOF[Np] = {"wTOF/expsigma/El", "wTOF/expsigma/Mu", "wTOF/expsigma/Pi",
                                                            "wTOF/expsigma/Ka", "wTOF/expsigma/Pr", "wTOF/expsigma/De",
                                                            "wTOF/expsigma/Tr", "wTOF/expsigma/He", "wTOF/expsigma/Al"};
  static constexpr std::string_view hnsigmaWithTOF[Np] = {"wTOF/nsigma/El", "wTOF/nsigma/Mu", "wTOF/nsigma/Pi",
                                                          "wTOF/nsigma/Ka", "wTOF/nsigma/Pr", "wTOF/nsigma/De",
                                                          "wTOF/nsigma/Tr", "wTOF/nsigma/He", "wTOF/nsigma/Al"};
  static constexpr std::string_view hnsigmaptWithTOF[Np] = {"wTOF/nsigmapt/El", "wTOF/nsigmapt/Mu", "wTOF/nsigmapt/Pi",
                                                            "wTOF/nsigmapt/Ka", "wTOF/nsigmapt/Pr", "wTOF/nsigmapt/De",
                                                            "wTOF/nsigmapt/Tr", "wTOF/nsigmapt/He", "wTOF/nsigmapt/Al"};
  static constexpr std::string_view hnsigmaposptWithTOF[Np] = {"wTOF/nsigmapospt/El", "wTOF/nsigmapospt/Mu", "wTOF/nsigmapospt/Pi",
                                                               "wTOF/nsigmapospt/Ka", "wTOF/nsigmapospt/Pr", "wTOF/nsigmapospt/De",
                                                               "wTOF/nsigmapospt/Tr", "wTOF/nsigmapospt/He", "wTOF/nsigmapospt/Al"};
  static constexpr std::string_view hnsigmanegptWithTOF[Np] = {"wTOF/nsigmanegpt/El", "wTOF/nsigmanegpt/Mu", "wTOF/nsigmanegpt/Pi",
                                                               "wTOF/nsigmanegpt/Ka", "wTOF/nsigmanegpt/Pr", "wTOF/nsigmanegpt/De",
                                                               "wTOF/nsigmanegpt/Tr", "wTOF/nsigmanegpt/He", "wTOF/nsigmanegpt/Al"};
  static constexpr std::string_view hsignalWithTOF[Np] = {"wTOF/signal/El", "wTOF/signal/Mu", "wTOF/signal/Pi",
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
    histos.add(hnsigmapt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmapospt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmanegpt[id].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});

    if (!enableFullHistos) { // Enabling only NSigma for tiny tables
      return;
    }

    // Exp signal
    const AxisSpec expAxis{1000, 0, 1000, Form("d#it{E}/d#it{x}_(%s) A.U.", pT[id])};
    histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)", pT[id])};
    histos.add(hexpected_diff[id].data(), "", kTH2F, {pAxis, deltaAxis});
    histos.add(hexpected_diffptpos[id].data(), "Positive", kTH2F, {pAxis, deltaAxis});
    histos.add(hexpected_diffptneg[id].data(), "Negative", kTH2F, {pAxis, deltaAxis});

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TPC}(%s)", pT[id])};
    histos.add(hexpsigma[id].data(), "", kTH2F, {pAxis, expSigmaAxis});

    if (!enableTOFHistos) { // Returning if the plots with TOF are not requested
      return;
    }

    histos.add(hexpectedWithTOF[id].data(), "With TOF", kTH2F, {pAxis, expAxis});
    histos.add(hexpected_diffWithTOF[id].data(), "With TOF", kTH2F, {pAxis, deltaAxis});
    histos.add(hexpected_diffptposWithTOF[id].data(), "With TOF Positive", kTH2F, {pAxis, deltaAxis});
    histos.add(hexpected_diffptnegWithTOF[id].data(), "With TOF Negative", kTH2F, {pAxis, deltaAxis});
    histos.add(hexpsigmaWithTOF[id].data(), "With TOF", kTH2F, {pAxis, expSigmaAxis});
    histos.add(hnsigmaWithTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigmaptWithTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmaposptWithTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmanegptWithTOF[id].data(), Form("With TOF %s", axisTitle), kTH2F, {ptAxis, nSigmaAxis});

    const AxisSpec dedxAxis{1000, 0, 1000, "d#it{E}/d#it{x} A.U."};
    histos.add(hsignalWithTOF[id].data(), "With TOF", kTH2F, {pAxis, dedxAxis});
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec multAxis{1000, 0.f, 1000.f, "Track multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec phiAxis{100, 0, TMath::TwoPi(), "#it{#phi}"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec pAxisPosNeg{nBinsP, -maxP, maxP, "Signed #it{p} (GeV/#it{c})"};
    const AxisSpec ptAxisPosNeg{nBinsP, -maxP, maxP, "Signed #it{p}_{T} (GeV/#it{c})"};
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
        if (t.sign() > 0) {
          histos.fill(HIST(hexpected_diffptpos[id]), t.pt(), diff);
        } else {
          histos.fill(HIST(hexpected_diffptneg[id]), t.pt(), diff);
        }
        histos.fill(HIST(hexpsigma[id]), t.tpcInnerParam(), o2::aod::pidutils::tpcExpSigma<id>(t));
        if constexpr (fillWithTOFHistograms) {
          if (std::abs(o2::aod::pidutils::tofNSigma<id>(t)) < 3.f) {
            histos.fill(HIST(hexpectedWithTOF[id]), t.tpcInnerParam(), t.tpcSignal() - diff);
            histos.fill(HIST(hexpected_diffWithTOF[id]), t.tpcInnerParam(), diff);
            histos.fill(HIST(hexpsigmaWithTOF[id]), t.p(), o2::aod::pidutils::tpcExpSigma<id>(t));
          }
        }
      }
      if constexpr (fillWithTOFHistograms) { // Filling nsigma (common to full and tiny)
        const auto& nsigmatof = o2::aod::pidutils::tofNSigma<id>(t);
        if (std::abs(nsigmatof) < 3.f) {
          histos.fill(HIST(hnsigmaWithTOF[id]), t.p(), nsigma);
          histos.fill(HIST(hnsigmaptWithTOF[id]), t.pt(), nsigma);
          histos.fill(HIST(hsignalWithTOF[id]), t.tpcInnerParam(), t.tpcSignal());
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
