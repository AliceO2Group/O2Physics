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
/// \file   qaPIDTOFBeta.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce the TOF QA plots for Beta
///

#include "Framework/AnalysisTask.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "qaPIDTOF.h"
#include "Framework/runDataProcessing.h"

struct tofPidBetaQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hdelta[Np] = {"delta/El", "delta/Mu", "delta/Pi",
                                                  "delta/Ka", "delta/Pr", "delta/De",
                                                  "delta/Tr", "delta/He", "delta/Al"};
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
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<bool> applyTrackCut{"applyTrackCut", false, "Flag to apply standard track cuts"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal"};
    const AxisSpec betaAxis{nBinsBeta, minBeta, maxBeta, "TOF #beta"};
    const AxisSpec massAxis{1000, 0, 3, "TOF mass (GeV/#it{c}^{2})"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    const AxisSpec pAxisPosNeg{2 * nBinsP, -maxP, maxP, "#it{p}/z (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
    }

    // Event properties
    histos.add("event/tofsignal", "", HistType::kTH2F, {pAxis, tofAxis});
    histos.add("event/tofmass", "TOF mass", HistType::kTH1F, {massAxis});
    histos.add("event/tofmassEvTimeTOF", "TOF mass Ev. Time TOF", HistType::kTH2F, {pAxis, massAxis});
    histos.add("event/tofmassEvTimeTOFOnly", "TOF mass Ev. Time TOF Only", HistType::kTH2F, {pAxis, massAxis});
    histos.add("event/tofmassEvTimeT0AC", "TOF mass Ev. Time T0AC", HistType::kTH2F, {pAxis, massAxis});
    histos.add("event/tofmassEvTimeT0ACOnly", "TOF mass Ev. Time T0AC Only", HistType::kTH2F, {pAxis, massAxis});
    histos.add("event/tofbeta", "", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/tofbetaEvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("event/signedtofbeta", "", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/signedtofbetaEvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("event/eta", "", HistType::kTH1F, {etaAxis});
    histos.add("event/length", "", HistType::kTH1F, {lAxis});
    histos.add("event/pt", "", HistType::kTH1F, {ptAxis});
    histos.add("event/p", "", HistType::kTH1F, {pAxis});
    auto h = histos.add<TH1>("event/evsel", "", kTH1F, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed mult.");
    h->GetXaxis()->SetBinLabel(4, "Passed vtx Z");

    h = histos.add<TH1>("event/trackselection", "", kTH1F, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "hasTOF");
    h->GetXaxis()->SetBinLabel(3, "isGlobalTrack");
  }

  template <uint8_t i, typename T>
  void fillParticleHistos(const T& t, const float tof, const float exp_diff, const float nsigma)
  {
    histos.fill(HIST(hexpected[i]), t.p(), tof - exp_diff);
    histos.fill(HIST(hdelta[i]), t.p(), exp_diff);
    histos.fill(HIST(hnsigma[i]), t.p(), nsigma);
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  void process(CollisionCandidate const& collision,
               soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFbeta, aod::pidTOFmass, aod::TrackSelection, aod::TOFSignal, aod::pidEvTimeFlags> const& tracks)
  {

    histos.fill(HIST("event/evsel"), 1);
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return;
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
      return;
    }

    histos.fill(HIST("event/evsel"), 4);

    for (auto const& track : tracks) {
      histos.fill(HIST("event/trackselection"), 1.f);
      if (!track.hasTOF()) { // Skipping tracks without TOF
        continue;
      }
      histos.fill(HIST("event/trackselection"), 2.f);
      if (applyTrackCut && !track.isGlobalTrack()) {
        continue;
      }
      histos.fill(HIST("event/trackselection"), 3.f);
      if (track.isEvTimeTOF()) {
        histos.fill(HIST("event/tofmassEvTimeTOF"), track.p(), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeTOF"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeTOF"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
        histos.fill(HIST("event/tofmassEvTimeTOFOnly"), track.p(), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeTOFOnly"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeTOFOnly"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC()) {
        histos.fill(HIST("event/tofmassEvTimeT0AC"), track.p(), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeT0AC"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeT0AC"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
        histos.fill(HIST("event/tofmassEvTimeT0ACOnly"), track.p(), track.mass());
        histos.fill(HIST("event/tofbetaEvTimeT0ACOnly"), track.p(), track.beta());
        histos.fill(HIST("event/signedtofbetaEvTimeT0ACOnly"), track.p() * track.sign(), track.beta());
      }
      histos.fill(HIST("event/tofmass"), track.p(), track.mass());
      histos.fill(HIST("event/tofbeta"), track.p(), track.beta());
      histos.fill(HIST("event/signedtofbeta"), track.p() * track.sign(), track.beta());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofPidBetaQa>(cfgc)};
}
