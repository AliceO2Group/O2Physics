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

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Task to produce the TOF Beta QA plots
struct tofPidBetaQa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> splitTrdTracks{"splitTrdTracks", false, "Flag to fill histograms for tracks with TRD match"};
  Configurable<bool> splitSignalPerCharge{"splitSignalPerCharge", true, "Split the signal per charge (reduces memory footprint if off)"};
  Configurable<bool> splitSignalPerEvTime{"splitSignalPerEvTime", true, "Split the signal per event time (reduces memory footprint if off)"};
  Configurable<int> lastTrdLayerForTrdMatch{"lastTrdLayerForTrdMatch", 5, "Last TRD layer to consider for TRD match"};

  ConfigurableAxis tofMassBins{"tofMassBins", {1000, 0, 3.f}, "Binning in the TOF mass plot"};
  ConfigurableAxis tofBetaBins{"tofBetaBins", {4000, 0, 2.f}, "Binning in the TOF beta plot"};
  ConfigurableAxis trackLengthBins{"trackLengthBins", {100, 0, 1000.f}, "Binning in track length plot"};
  Configurable<bool> requireGoodMatchTracks{"requireGoodMatchTracks", false, "Require good match tracks"};
  Configurable<float> mMaxTOFChi2{"maxTOFChi2", 3.f, "Maximum TOF Chi2"};
  Configurable<float> mEtaWindow{"etaWindow", 0.8f, "Window in eta for tracks"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal"};
    const AxisSpec betaAxis{tofBetaBins, "TOF #beta"};
    const AxisSpec massAxis{tofMassBins, "TOF mass (GeV/#it{c}^{2})"};
    const AxisSpec trdAxis{10, -0.5, 9.5, "Last TRD cluster"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec lAxis{trackLengthBins, "Track length (cm)"};
    const AxisSpec tofChi2Axis{1000, 0, 20, "TOF residual (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    const AxisSpec pAxisPosNeg{2 * nBinsP, -maxP, maxP, "signed #it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
    }

    // Event properties
    histos.add("event/tofsignal", "", HistType::kTH2F, {pAxis, tofAxis});
    const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};

    // TOF mass
    if (splitSignalPerCharge) {
      histos.add("tofmass/inclusive", "", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
      if (splitSignalPerEvTime) {
        histos.add("tofmass/EvTimeTOF", "Ev. Time TOF", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
        histos.add("tofmass/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
        histos.add("tofmass/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
        histos.add("tofmass/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
      }
      if (splitTrdTracks) {
        histos.add("tofmass/trd/inclusive", "(hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofmass/trd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
          histos.add("tofmass/trd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
          histos.add("tofmass/trd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
          histos.add("tofmass/trd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
        }
        histos.add("tofmass/notrd/inclusive", "(hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofmass/notrd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
          histos.add("tofmass/notrd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
          histos.add("tofmass/notrd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
          histos.add("tofmass/notrd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH3F, {pAxis, massAxis, chargeAxis});
        }
      }
    } else {
      histos.add("tofmass/inclusive", "", HistType::kTH2F, {pAxis, massAxis});
      if (splitSignalPerEvTime) {
        histos.add("tofmass/EvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxis, massAxis});
        histos.add("tofmass/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxis, massAxis});
        histos.add("tofmass/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxis, massAxis});
        histos.add("tofmass/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxis, massAxis});
      }
      if (splitTrdTracks) {
        histos.add("tofmass/trd/inclusive", "(hasTRD)", HistType::kTH2F, {pAxis, massAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofmass/trd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
          histos.add("tofmass/trd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
          histos.add("tofmass/trd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
          histos.add("tofmass/trd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
        }
        histos.add("tofmass/notrd/inclusive", "(hasTRD)", HistType::kTH2F, {pAxis, massAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofmass/notrd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
          histos.add("tofmass/notrd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
          histos.add("tofmass/notrd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
          histos.add("tofmass/notrd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH2F, {pAxis, massAxis});
        }
      }
    }

    // TOF beta
    if (splitSignalPerCharge) {
      histos.add("tofbeta/inclusive", "", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
      if (splitSignalPerEvTime) {
        histos.add("tofbeta/EvTimeTOF", "Ev. Time TOF", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
        histos.add("tofbeta/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
        histos.add("tofbeta/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
        histos.add("tofbeta/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
      }
      if (splitTrdTracks) {
        histos.add("tofbeta/trd/inclusive", "(hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofbeta/trd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
          histos.add("tofbeta/trd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
          histos.add("tofbeta/trd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
          histos.add("tofbeta/trd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
        }
        histos.add("tofbeta/notrd/inclusive", "(hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofbeta/notrd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
          histos.add("tofbeta/notrd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
          histos.add("tofbeta/notrd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
          histos.add("tofbeta/notrd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH3F, {pAxis, betaAxis, chargeAxis});
        }
      }
    } else {
      histos.add("tofbeta/inclusive", "", HistType::kTH2F, {pAxisPosNeg, betaAxis});
      if (splitSignalPerEvTime) {
        histos.add("tofbeta/EvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxisPosNeg, betaAxis});
        histos.add("tofbeta/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
        histos.add("tofbeta/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxisPosNeg, betaAxis});
        histos.add("tofbeta/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxisPosNeg, betaAxis});
      }
      if (splitTrdTracks) {
        histos.add("tofbeta/trd/inclusive", "(hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofbeta/trd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
          histos.add("tofbeta/trd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
          histos.add("tofbeta/trd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
          histos.add("tofbeta/trd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
        }
        histos.add("tofbeta/notrd/inclusive", "(hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
        if (splitSignalPerEvTime) {
          histos.add("tofbeta/notrd/EvTimeTOF", "Ev. Time TOF (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
          histos.add("tofbeta/notrd/EvTimeTOFOnly", "Ev. Time TOF Only (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
          histos.add("tofbeta/notrd/EvTimeT0AC", "Ev. Time T0AC (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
          histos.add("tofbeta/notrd/EvTimeT0ACOnly", "Ev. Time T0AC Only (hasTRD)", HistType::kTH2F, {pAxisPosNeg, betaAxis});
        }
      }
    }

    histos.add("event/tofchi2", "", HistType::kTH1F, {tofChi2Axis});
    histos.add("event/eta", "", HistType::kTH1F, {etaAxis});
    histos.add("event/length", "", HistType::kTH1F, {lAxis});
    if (splitTrdTracks) {
      histos.add("event/trd/length", "", HistType::kTH2F, {lAxis, trdAxis});
      histos.add("event/notrd/length", "", HistType::kTH1F, {lAxis});
    }
    histos.add("event/pt", "", HistType::kTH1F, {ptAxis});
    histos.add("event/p", "", HistType::kTH1F, {pAxis});
    auto h = histos.add<TH1>("event/evsel", "", kTH1F, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed vtx Z");

    h = histos.add<TH1>("event/trackselection", "", kTH1F, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "hasTOF");
    h->GetXaxis()->SetBinLabel(3, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(4, TString::Format("TOF chi2 < %.2f", mMaxTOFChi2.value));
  }

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (o2::aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (o2::aod::evsel::sel8 == true));
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireInAcceptanceTracksInFilter());
  Filter etaFilter = (nabs(o2::aod::track::eta) < mEtaWindow);

  using CollisionCandidate = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                    aod::pidTOFbeta, aod::pidTOFmass,
                                    aod::pidEvTimeFlags, aod::TOFSignal, aod::TOFEvTime,
                                    aod::pidTOFFlags>;
  void process(CollisionCandidate const& collision,
               soa::Filtered<TrackCandidates> const& tracks)
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

    if (std::abs(collision.posZ()) > 10.f) {
      return;
    }

    histos.fill(HIST("event/evsel"), 3);

    for (auto const& track : tracks) {
      histos.fill(HIST("event/trackselection"), 1.f);
      if (!track.hasTOF()) { // Skipping tracks without TOF
        continue;
      }
      histos.fill(HIST("event/trackselection"), 2.f);
      if (!track.isGlobalTrack()) {
        continue;
      }
      histos.fill(HIST("event/trackselection"), 3.f);
      if (track.tofChi2() > mMaxTOFChi2) { // Skipping tracks with large Chi2
        continue;
      }
      histos.fill(HIST("event/trackselection"), 4.f);

      if (splitSignalPerCharge) {
        histos.fill(HIST("tofmass/inclusive"), track.p(), track.mass(), track.sign());
        histos.fill(HIST("tofbeta/inclusive"), track.p(), track.beta(), track.sign());
        if (splitSignalPerEvTime) {
          if (track.isEvTimeTOF()) {
            histos.fill(HIST("tofmass/EvTimeTOF"), track.p(), track.mass(), track.sign());
            histos.fill(HIST("tofbeta/EvTimeTOF"), track.p(), track.beta(), track.sign());
          }
          if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
            histos.fill(HIST("tofmass/EvTimeTOFOnly"), track.p(), track.mass(), track.sign());
            histos.fill(HIST("tofbeta/EvTimeTOFOnly"), track.p(), track.beta(), track.sign());
          }
          if (track.isEvTimeT0AC()) {
            histos.fill(HIST("tofmass/EvTimeT0AC"), track.p(), track.mass(), track.sign());
            histos.fill(HIST("tofbeta/EvTimeT0AC"), track.p(), track.beta(), track.sign());
          }
          if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
            histos.fill(HIST("tofmass/EvTimeT0ACOnly"), track.p(), track.mass(), track.sign());
            histos.fill(HIST("tofbeta/EvTimeT0ACOnly"), track.p(), track.beta(), track.sign());
          }
        }
      } else {
        histos.fill(HIST("tofmass/inclusive"), track.p(), track.mass());
        histos.fill(HIST("tofbeta/inclusive"), track.p(), track.beta());
        if (splitSignalPerEvTime) {
          if (track.isEvTimeTOF()) {
            histos.fill(HIST("tofmass/EvTimeTOF"), track.p(), track.mass());
            histos.fill(HIST("tofbeta/EvTimeTOF"), track.p(), track.beta());
          }
          if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
            histos.fill(HIST("tofmass/EvTimeTOFOnly"), track.p(), track.mass());
            histos.fill(HIST("tofbeta/EvTimeTOFOnly"), track.p(), track.beta());
          }
          if (track.isEvTimeT0AC()) {
            histos.fill(HIST("tofmass/EvTimeT0AC"), track.p(), track.mass());
            histos.fill(HIST("tofbeta/EvTimeT0AC"), track.p(), track.beta());
          }
          if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
            histos.fill(HIST("tofmass/EvTimeT0ACOnly"), track.p(), track.mass());
            histos.fill(HIST("tofbeta/EvTimeT0ACOnly"), track.p(), track.beta());
          }
        }
      }
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/tofchi2"), track.tofChi2());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/tofsignal"), track.p(), track.tofSignal());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());

      if (!splitTrdTracks) { // If splitting of TRD tracks is not enabled, skip
        continue;
      }

      if (!track.hasTRD()) {
        histos.fill(HIST("event/notrd/length"), track.length());
        if (splitSignalPerCharge) {
          histos.fill(HIST("tofmass/notrd/inclusive"), track.p(), track.mass(), track.sign());
          histos.fill(HIST("tofbeta/notrd/inclusive"), track.p(), track.beta(), track.sign());
          if (splitSignalPerEvTime) {
            if (track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/notrd/EvTimeTOF"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/notrd/EvTimeTOF"), track.p(), track.beta(), track.sign());
            }
            if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/notrd/EvTimeTOFOnly"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/notrd/EvTimeTOFOnly"), track.p(), track.beta(), track.sign());
            }
            if (track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/notrd/EvTimeT0AC"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/notrd/EvTimeT0AC"), track.p(), track.beta(), track.sign());
            }
            if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/notrd/EvTimeT0ACOnly"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/notrd/EvTimeT0ACOnly"), track.p(), track.beta(), track.sign());
            }
          }
        } else {
          const float signedp = track.p() * track.sign();
          histos.fill(HIST("tofmass/notrd/inclusive"), signedp, track.mass());
          histos.fill(HIST("tofbeta/notrd/inclusive"), signedp, track.beta());
          if (splitSignalPerEvTime) {
            if (track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/notrd/EvTimeTOF"), signedp, track.mass());
              histos.fill(HIST("tofbeta/notrd/EvTimeTOF"), signedp, track.beta());
            }
            if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/notrd/EvTimeTOFOnly"), signedp, track.mass());
              histos.fill(HIST("tofbeta/notrd/EvTimeTOFOnly"), signedp, track.beta());
            }
            if (track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/notrd/EvTimeT0AC"), signedp, track.mass());
              histos.fill(HIST("tofbeta/notrd/EvTimeT0AC"), signedp, track.beta());
            }
            if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/notrd/EvTimeT0ACOnly"), signedp, track.mass());
              histos.fill(HIST("tofbeta/notrd/EvTimeT0ACOnly"), signedp, track.beta());
            }
          }
        }
      } else {

        int lastLayer = 0;
        for (int l = 7; l >= 0; l--) {
          if (track.trdPattern() & (1 << l)) {
            lastLayer = l;
            break;
          }
        }

        histos.fill(HIST("event/trd/length"), track.length(), lastLayer);
        if (lastLayer < lastTrdLayerForTrdMatch) {
          continue;
        }
        if (splitSignalPerCharge) {
          histos.fill(HIST("tofmass/trd/inclusive"), track.p(), track.mass(), track.sign());
          histos.fill(HIST("tofbeta/trd/inclusive"), track.p(), track.beta(), track.sign());
          if (splitSignalPerEvTime) {
            if (track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/trd/EvTimeTOF"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/trd/EvTimeTOF"), track.p(), track.beta(), track.sign());
            }
            if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/trd/EvTimeTOFOnly"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/trd/EvTimeTOFOnly"), track.p(), track.beta(), track.sign());
            }
            if (track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/trd/EvTimeT0AC"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/trd/EvTimeT0AC"), track.p(), track.beta(), track.sign());
            }
            if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/trd/EvTimeT0ACOnly"), track.p(), track.mass(), track.sign());
              histos.fill(HIST("tofbeta/trd/EvTimeT0ACOnly"), track.p(), track.beta(), track.sign());
            }
          }
        } else {
          const float signedp = track.p() * track.sign();
          histos.fill(HIST("tofmass/trd/inclusive"), signedp, track.mass());
          histos.fill(HIST("tofbeta/trd/inclusive"), signedp, track.beta());
          if (splitSignalPerEvTime) {
            if (track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/trd/EvTimeTOF"), signedp, track.mass());
              histos.fill(HIST("tofbeta/trd/EvTimeTOF"), signedp, track.beta());
            }
            if (track.isEvTimeTOF() && !track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/trd/EvTimeTOFOnly"), signedp, track.mass());
              histos.fill(HIST("tofbeta/trd/EvTimeTOFOnly"), signedp, track.beta());
            }
            if (track.isEvTimeT0AC()) {
              histos.fill(HIST("tofmass/trd/EvTimeT0AC"), signedp, track.mass());
              histos.fill(HIST("tofbeta/trd/EvTimeT0AC"), signedp, track.beta());
            }
            if (track.isEvTimeT0AC() && !track.isEvTimeTOF()) {
              histos.fill(HIST("tofmass/trd/EvTimeT0ACOnly"), signedp, track.mass());
              histos.fill(HIST("tofbeta/trd/EvTimeT0ACOnly"), signedp, track.beta());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tofPidBetaQa>(cfgc)}; }
