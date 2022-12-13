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
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/FT0Corrected.h"

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
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<bool> applyTrackCut{"applyTrackCut", true, "Flag to apply standard track cuts"};
  Configurable<bool> splitTrdTracks{"splitTrdTracks", false, "Flag to fill histograms for tracks with TRD match"};
  ConfigurableAxis tofMassBins{"tofMassBins", {1000, 0, 3.f}, "Binning in the TOF mass plot"};
  ConfigurableAxis tofBetaBins{"tofBetaBins", {4000, 0, 2.f}, "Binning in the TOF beta plot"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal"};
    const AxisSpec betaAxis{tofBetaBins, "TOF #beta"};
    const AxisSpec massAxis{tofMassBins, "TOF mass (GeV/#it{c}^{2})"};
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
    histos.add("tofmass/inclusive", "", HistType::kTH2F, {pAxis, massAxis});
    histos.add("tofmass/EvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxis, massAxis});
    histos.add("tofmass/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxis, massAxis});
    histos.add("tofmass/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxis, massAxis});
    histos.add("tofmass/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxis, massAxis});
    if (splitTrdTracks) {
      histos.add("tofmass/trd/inclusive", "(hasTRD)", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/trd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/trd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/trd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/trd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
    }

    histos.add("tofbeta/inclusive", "", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("tofbeta/EvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("tofbeta/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("tofbeta/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxis, betaAxis});
    histos.add("tofbeta/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxis, betaAxis});
    if (splitTrdTracks) {
      histos.add("tofbeta/trd/inclusive", "(hasTRD)", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/trd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/trd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/trd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/trd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH2F, {pAxis, betaAxis});
    }

    histos.add("signedtofbeta/inclusive", "", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("signedtofbeta/EvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("signedtofbeta/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("signedtofbeta/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    histos.add("signedtofbeta/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    if (splitTrdTracks) {
      histos.add("signedtofbeta/trd/inclusive", " (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
      histos.add("signedtofbeta/trd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
      histos.add("signedtofbeta/trd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
      histos.add("signedtofbeta/trd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
      histos.add("signedtofbeta/trd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
    }

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
        histos.fill(HIST("tofmass/EvTimeTOF"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/EvTimeTOF"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/EvTimeTOF"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
        histos.fill(HIST("tofmass/EvTimeTOFOnly"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/EvTimeTOFOnly"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/EvTimeTOFOnly"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC()) {
        histos.fill(HIST("tofmass/EvTimeT0AC"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/EvTimeT0AC"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/EvTimeT0AC"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
        histos.fill(HIST("tofmass/EvTimeT0ACOnly"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/EvTimeT0ACOnly"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/EvTimeT0ACOnly"), track.p() * track.sign(), track.beta());
      }
      histos.fill(HIST("tofmass/inclusive"), track.p(), track.mass());
      histos.fill(HIST("tofbeta/inclusive"), track.p(), track.beta());
      histos.fill(HIST("signedtofbeta/inclusive"), track.p() * track.sign(), track.beta());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
      if (!splitTrdTracks || !track.hasTRD()) {
        continue;
      }
      if (track.isEvTimeTOF()) {
        histos.fill(HIST("tofmass/trd/EvTimeTOF"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/trd/EvTimeTOF"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/trd/EvTimeTOF"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
        histos.fill(HIST("tofmass/trd/EvTimeTOFOnly"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/trd/EvTimeTOFOnly"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/trd/EvTimeTOFOnly"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC()) {
        histos.fill(HIST("tofmass/trd/EvTimeT0AC"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/trd/EvTimeT0AC"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/trd/EvTimeT0AC"), track.p() * track.sign(), track.beta());
      }
      if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
        histos.fill(HIST("tofmass/trd/EvTimeT0ACOnly"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/trd/EvTimeT0ACOnly"), track.p(), track.beta());
        histos.fill(HIST("signedtofbeta/trd/EvTimeT0ACOnly"), track.p() * track.sign(), track.beta());
      }
      histos.fill(HIST("tofmass/trd/inclusive"), track.p(), track.mass());
      histos.fill(HIST("tofbeta/trd/inclusive"), track.p(), track.beta());
      histos.fill(HIST("signedtofbeta/trd/inclusive"), track.p() * track.sign(), track.beta());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofPidBetaQa>(cfgc)};
}
