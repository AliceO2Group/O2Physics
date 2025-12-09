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
/// \file   tofOfflineCalib.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  14/02/2023
/// \brief  Task to produce calibration objects for the TOF. Based on AO2D or TOF skimmed data
///

#include "tofSkimsTableCreator.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

struct tofOfflineCalib {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::TOFEvTime, aod::EvTimeTOFOnly, aod::TOFSignal, aod::pidEvTimeFlags,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                         aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::FT0sCorrected>;

  // Tables to be produced
  Produces<o2::aod::DeltaTOF> tableRow;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> makeTable{"makeTable", false, "Make an output table"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 0.1, "Fractions of events to keep"};
  Configurable<float> pRefMin{"pRefMin", 0.6, "Reference momentum minimum"};
  Configurable<float> pRefMax{"pRefMax", 0.7, "Reference momentum maximum"};
  Configurable<float> maxTOFChi2{"maxTOFChi2", 5.f, "Maximum TOF Chi2 to accept tracks"};
  Configurable<float> deltatTh{"deltatTh", 500, "Threshold in DeltaT to accept reference tracks"};

  std::shared_ptr<TH2> deltaVsP;
  std::shared_ptr<TH2> deltaVsPHighChi2;
  std::shared_ptr<TH1> hGood;
  std::shared_ptr<TH1> hBad;
  std::shared_ptr<TH1> hGoodRefWithTRD;
  std::shared_ptr<TH1> hBadRefWithTRD;

  unsigned int randomSeed = 0;
  void init(o2::framework::InitContext&)
  {
    randomSeed = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    auto h = histos.add<TH1>("events", "Events", kTH1D, {{10, 0, 10, "Event selection"}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Events sampled");

    switch (applyEvSel.value) {
      case 0:
      case 1:
      case 2:
        break;
      default:
        LOG(fatal) << "Invalid event selection flag: " << applyEvSel.value;
        break;
    }
    switch (trackSelection.value) {
      case 0:
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
        break;
      default:
        LOG(fatal) << "Invalid track selection flag: " << trackSelection.value;
        break;
    }
  }

  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<soa::Filtered<Trks>::iterator, pid>;

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (o2::aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (o2::aod::evsel::sel8 == true));
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireInAcceptanceTracksInFilter());

  int lastRun = -1;
  void process(soa::Filtered<Coll>::iterator const& collision,
               soa::Filtered<Trks> const& tracks,
               aod::BCs const&)
  {
    histos.fill(HIST("events"), 0.5);
    if (fractionOfEvents < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      return;
    }
    histos.fill(HIST("events"), 1.5);

    if (makeTable) { // Reserve room for the output table
      tableRow.reserve(tracks.size());
    }

    // Handle event times
    float evTimeT0AC = 0.f;
    float evTimeT0ACErr = 0.f;
    if (collision.t0ACValid()) {
      evTimeT0AC = collision.t0AC() * 1000.f;
      evTimeT0ACErr = collision.t0resolution() * 1000.f;
    }
    float evTimeT0A = 0.f;
    if (collision.t0ACorrectedValid()) {
      evTimeT0A = collision.t0ACorrected() * 1000.f;
    }
    float evTimeT0C = 0.f;
    if (collision.t0CCorrectedValid()) {
      evTimeT0C = collision.t0CCorrected() * 1000.f;
    }

    if (lastRun != collision.bc().runNumber()) {
      lastRun = collision.bc().runNumber();
      const AxisSpec doubleDeltaAxis{1500, -3000, 3000, "#Deltat_{#pi} - #Deltat_{#pi}^{ref} (ps)"};
      const AxisSpec pTAxis{200, 0, 5, "#it{p}_{T} (GeV/#it{c})"};
      hGood = histos.add<TH1>(Form("Run%i/hGood", lastRun), "Good", kTH1D, {doubleDeltaAxis});
      hBad = histos.add<TH1>(Form("Run%i/hBad", lastRun), "Bad", kTH1D, {doubleDeltaAxis});
      hGoodRefWithTRD = histos.add<TH1>(Form("Run%i/hGoodRefWithTRD", lastRun), "Good", kTH1D, {doubleDeltaAxis});
      hBadRefWithTRD = histos.add<TH1>(Form("Run%i/hBadRefWithTRD", lastRun), "Bad", kTH1D, {doubleDeltaAxis});
      deltaVsP = histos.add<TH2>(Form("Run%i/deltaVsP", lastRun), "Low Chi2", kTH2F, {pTAxis, doubleDeltaAxis});
      deltaVsPHighChi2 = histos.add<TH2>(Form("Run%i/deltaVsPHighChi2", lastRun), "High Chi2", kTH2F, {pTAxis, doubleDeltaAxis});
    }

    int8_t lastTRDLayer = -1;
    for (auto& track1 : tracks) {
      if (!track1.hasTOF()) {
        continue;
      }
      // Selecting good reference
      const float& texp1Pi = track1.tofExpSignalPi(track1.tofSignal() - track1.tofEvTime());
      const float& delta1Pi = track1.tofSignal() - texp1Pi;
      if (track1.p() < pRefMin || track1.p() > pRefMax || track1.tofChi2() > maxTOFChi2 || fabs(delta1Pi) > deltatTh) {
        continue;
      }
      for (auto& track2 : tracks) {
        if (!track2.hasTOF()) {
          continue;
        }
        if (track1.globalIndex() == track2.globalIndex()) { // Skipping the same track
          continue;
        }
        const float& texp2El = track2.tofExpSignalEl(track2.tofSignal() - track2.tofEvTime());
        const float& texp2Mu = track2.tofExpSignalMu(track2.tofSignal() - track2.tofEvTime());
        const float& texp2Pi = track2.tofExpSignalPi(track2.tofSignal() - track2.tofEvTime());
        const float& texp2Ka = track2.tofExpSignalKa(track2.tofSignal() - track2.tofEvTime());
        const float& texp2Pr = track2.tofExpSignalPr(track2.tofSignal() - track2.tofEvTime());
        const float& delta2Pi = track2.tofSignal() - texp2Pi;
        if (track2.tofChi2() < maxTOFChi2) {
          deltaVsP->Fill(track2.p(), delta2Pi - delta1Pi);
        } else if (track2.tofChi2() > maxTOFChi2) {
          deltaVsPHighChi2->Fill(track2.p(), delta2Pi - delta1Pi);
        }
        if (track2.p() > pRefMin && track2.p() < pRefMax) {
          if (track1.hasTRD()) {
            if (track2.tofChi2() < maxTOFChi2) {
              hGoodRefWithTRD->Fill(delta2Pi - delta1Pi);
            } else if (track2.tofChi2() > maxTOFChi2 + 2) {
              hBadRefWithTRD->Fill(delta2Pi - delta1Pi);
            }
          } else {
            if (track2.tofChi2() < maxTOFChi2) {
              hGood->Fill(delta2Pi - delta1Pi);
            } else if (track2.tofChi2() > maxTOFChi2 + 2) {
              hBad->Fill(delta2Pi - delta1Pi);
            }
          }
        }
        lastTRDLayer = -1;
        if (track2.hasTRD()) {
          for (int8_t l = 7; l >= 0; l--) {
            if (track2.trdPattern() & (1 << l)) {
              lastTRDLayer = l;
              break;
            }
          }
        }

        if (!makeTable) {
          continue;
        }
        tableRow(track2.collisionId(),
                 track2.p(),
                 track1.p() - track2.p(),
                 track2.pt() * track2.sign(),
                 track1.pt() - track2.pt(),
                 track2.eta(),
                 track1.eta() - track2.eta(),
                 track2.phi(),
                 track1.phi() - track2.phi(),
                 track2.tofSignal() - texp2El,
                 track2.tofSignal() - texp2Mu,
                 delta2Pi,
                 track2.tofSignal() - texp2Ka,
                 track2.tofSignal() - texp2Pr,
                 delta2Pi - delta1Pi,
                 track1.sign(),
                 track2.length(),
                 track2.tofChi2(),
                 track2.tpcSignal(),
                 track2.tofSignal(),
                 track2.evTimeTOF(),
                 track2.evTimeTOFErr(),
                 track2.evTimeTOFMult(),
                 evTimeT0A,
                 evTimeT0C,
                 evTimeT0AC,
                 evTimeT0ACErr,
                 collision.collisionTime(),
                 collision.collisionTimeRes(),
                 track2.tofFlags(),
                 lastTRDLayer);

        // float doubleDelta = delta2Pi - delta1Pi;

        // if (fabs(doubleDelta) < 3000) // Good electrons
        //   tout->Fill();
        // if (fabs(doubleDelta) < 3000) // Good pions
        //   tout->Fill();
        // if (fabs(doubleDelta) < 3000) // Good kaons
        //   tout->Fill();
        // if (fabs(doubleDelta) < 3000) // Good protons
        //   tout->Fill();
        // if (fabs(doubleDelta) < 3000) // Good deuterons
        //   tout->Fill();
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tofOfflineCalib>(cfgc)}; }
