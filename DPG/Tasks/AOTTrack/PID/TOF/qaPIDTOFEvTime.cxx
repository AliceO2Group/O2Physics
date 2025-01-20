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
/// \file   qaPIDTOFEvTime.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Tasks of the TOF PID quantities for the event times
///

#include "TEfficiency.h"
#include "THashList.h"

#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Framework/runDataProcessing.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFT0/Digit.h"

using namespace o2;
using namespace o2::framework;

struct tofPidCollisionTimeQa {
  ConfigurableAxis evTimeBins{"evTimeBins", {1000, -1000.f, 1000.f}, "Binning for the event time"};
  ConfigurableAxis evTimeDeltaBins{"evTimeDeltaBins", {1000, -1000.f, 1000.f}, "Binning for the delta between event times"};
  ConfigurableAxis evTimeResoBins{"evTimeResoBins", {1000, 0.f, 1000.f}, "Binning for the event time resolution"};
  ConfigurableAxis tofSignalBins{"tofSignalBins", {5000, 0.f, 100000.f}, "Binning for the TOF signal"};
  ConfigurableAxis pBins{"pBins", {200, 0.1f, 5.f}, "Binning for the momentum"};

  Configurable<int> nBinsMultiplicity{"nBinsMultiplicity", 1000, "Number of bins for the multiplicity"};
  Configurable<float> rangeMultiplicity{"rangeMultiplicity", 1000.f, "Range for the multiplicity"};
  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<float> minPReso{"minPReso", 1.4f, "Minimum momentum in range for the resolution plot"};
  Configurable<float> maxPReso{"maxPReso", 1.5f, "Maximum momentum in range for the resolution plot"};
  Configurable<bool> enableDebug{"enableDebug", false, "Add debug plots"};

  OutputObj<THashList> listEfficiency{"Efficiency"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    const AxisSpec evTimeAxis{evTimeBins, "Event time (ps)"};
    const AxisSpec evTimeDeltaAxis{evTimeDeltaBins, "Delta event time (ps)"};
    const AxisSpec multAxis{nBinsMultiplicity, 0, rangeMultiplicity, "Track multiplicity for TOF event time"};
    const AxisSpec evTimeResoAxis{evTimeResoBins, "Event time resolution (ps)"};
    const AxisSpec tofSignalAxis{tofSignalBins, "TOF signal (ps)"};
    AxisSpec pAxis{pBins, "#it{p} GeV/#it{c}"};
    AxisSpec ptAxis{pBins, "#it{p}_{T} GeV/#it{c}"};
    if (logAxis) {
      pAxis.makeLogarithmic();
      ptAxis.makeLogarithmic();
    }
    const AxisSpec collisionAxis{6000, -0.5f, 6000.f - .5f, "Collision index % 6000"};
    const AxisSpec massAxis{1000, 0, 3, "TOF mass (GeV/#it{c}^{2})"};
    const AxisSpec betaAxis{1000, 0, 1.5, "TOF #beta"};
    const AxisSpec deltaAxis{1000, -10000, 10000, "t-t_{ev}-t_{exp}(#pi) (ps)"};
    const AxisSpec lengthAxis{1000, 0, 600, "Track length (cm)"};

    if (doprocessData) {
      auto h = histos.add<TH1>("eventSelection", "eventSelection", kTH1F, {{10, 0, 10, "Cut passed"}});
      h->GetXaxis()->SetBinLabel(1, "Events read");
      h->GetXaxis()->SetBinLabel(2, "Event selection");
      h->GetXaxis()->SetBinLabel(3, "#sigma_{Ev. time} < 200 ps");
      h->GetXaxis()->SetBinLabel(4, "#sigma_{Ev. time} > 200 ps");
      h = histos.add<TH1>("trackSelection", "trackSelection", kTH1F, {{10, 0, 10, "Cut passed"}});
      h->GetXaxis()->SetBinLabel(1, "Tracks read");
      h->GetXaxis()->SetBinLabel(2, "Track selection");
      h->GetXaxis()->SetBinLabel(3, "hasITS");
      h->GetXaxis()->SetBinLabel(4, "hasTPC");
      h->GetXaxis()->SetBinLabel(5, "hasTOF");
      histos.add<TH1>("deltaEvTimeTOFT0A", "deltaEvTimeTOFT0A", kTH1F, {evTimeDeltaAxis})->GetXaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0A} (ps)");
      histos.add<TH1>("deltaEvTimeTOFT0C", "deltaEvTimeTOFT0C", kTH1F, {evTimeDeltaAxis})->GetXaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0C} (ps)");
      histos.add<TH1>("deltaEvTimeTOFT0AC", "deltaEvTimeTOFT0AC", kTH1F, {evTimeDeltaAxis})->GetXaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0AC} (ps)");
      auto h2 = histos.add<TH2>("deltaEvTimeTOFT0AvsT0C", "deltaEvTimeTOFT0AvsT0C", kTH2F, {evTimeDeltaAxis, evTimeDeltaAxis});
      h2->GetXaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0A} (ps)");
      h2->GetYaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0C} (ps)");
      h2 = histos.add<TH2>("deltaEvTimeTOFT0AvsTOF", "deltaEvTimeTOFT0AvsTOF", kTH2F, {evTimeAxis, evTimeDeltaAxis});
      h2->GetXaxis()->SetTitle("Ev. time_{TOF} (ps)");
      h2->GetYaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0A} (ps)");
      h2 = histos.add<TH2>("deltaEvTimeTOFT0CvsTOF", "deltaEvTimeTOFT0CvsTOF", kTH2F, {evTimeAxis, evTimeDeltaAxis});
      h2->GetXaxis()->SetTitle("Ev. time_{TOF} (ps)");
      h2->GetYaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0C} (ps)");
      h2 = histos.add<TH2>("deltaEvTimeTOFT0ACvsTOF", "deltaEvTimeTOFT0ACvsTOF", kTH2F, {evTimeAxis, evTimeDeltaAxis});
      h2->GetXaxis()->SetTitle("Ev. time_{TOF} (ps)");
      h2->GetYaxis()->SetTitle("Ev. time_{TOF} - Ev. time_{T0AC} (ps)");

      histos.add("eventTime", "eventTime", kTH1F, {evTimeAxis});
      histos.add("eventTimeReso", "eventTimeReso", kTH1F, {evTimeResoAxis});
      histos.add("eventTimeVsMult", "eventTimeVsMult", kTH2F, {multAxis, evTimeAxis});
      histos.add("eventTimeResoVsMult", "eventTimeResoVsMult", kTH2F, {multAxis, evTimeResoAxis});

      histos.add("eventTimeTOFMult", "eventTimeTOFMult", kTH1F, {multAxis});
      histos.add<TH1>("eventTimeTOF", "eventTimeTOF", kTH1F, {evTimeAxis})->GetXaxis()->SetTitle("Ev. time_{TOF} (ps)");
      histos.add<TH1>("eventTimeTOFReso", "eventTimeTOFReso", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("Ev. time_{TOF} resolution (ps)");
      histos.add<TH2>("eventTimeTOFVsMult", "eventTimeTOFVsMult", kTH2F, {multAxis, evTimeAxis})->GetYaxis()->SetTitle("Ev. time_{TOF} (ps)");
      histos.add<TH2>("eventTimeTOFResoVsMult", "eventTimeTOFResoVsMult", kTH2F, {multAxis, evTimeResoAxis})->GetYaxis()->SetTitle("Ev. time_{TOF} resolution (ps)");

      histos.add<TH1>("eventTimeT0A", "eventTimeT0A", kTH1F, {evTimeAxis})->GetXaxis()->SetTitle("T0A event time (ps)");
      histos.add<TH1>("eventTimeT0C", "eventTimeT0C", kTH1F, {evTimeAxis})->GetXaxis()->SetTitle("T0C event time (ps)");
      histos.add<TH1>("eventTimeT0AC", "eventTimeT0AC", kTH1F, {evTimeAxis})->GetXaxis()->SetTitle("T0AC event time (ps)");
      histos.add<TH1>("eventTimeT0ACReso", "eventTimeT0ACReso", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("T0AC event time resolution (ps)");

      histos.add<TH1>("collisionTime", "collisionTime", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("Collision time (ps)");
      histos.add<TH1>("collisionTimeRes", "collisionTimeRes", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("Collision time resolution (ps)");

      histos.add("tracks/p", "p", kTH1F, {pAxis});
      histos.add("tracks/pt", "pt", kTH1F, {ptAxis});
      histos.add("tracks/length", "length", kTH1F, {lengthAxis});

      histos.add("deltaVsMult/pi", Form("pi %.2f < #it{p} < %.2f", minPReso.value, maxPReso.value), kTH2F, {multAxis, deltaAxis});
      histos.add("deltaVsReso/pi", Form("pi %.2f < #it{p} < %.2f", minPReso.value, maxPReso.value), kTH2F, {evTimeResoAxis, deltaAxis});

      histos.add("tofbeta/inclusive", "", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/EvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/EvTimeT0AOnly", "Ev. Time T0A Only", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/EvTimeT0COnly", "Ev. Time T0A Only", HistType::kTH2F, {pAxis, betaAxis});
      histos.add("tofbeta/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxis, betaAxis});

      histos.add("tofmass/inclusive", "", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/EvTimeTOF", "Ev. Time TOF", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/EvTimeTOFOnly", "Ev. Time TOF Only", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/EvTimeT0AC", "Ev. Time T0AC", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/EvTimeT0AOnly", "Ev. Time T0A Only", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/EvTimeT0COnly", "Ev. Time T0C Only", HistType::kTH2F, {pAxis, massAxis});
      histos.add("tofmass/EvTimeT0ACOnly", "Ev. Time T0AC Only", HistType::kTH2F, {pAxis, massAxis});
      if (enableDebug) {
        histos.add("withtof/p", "p", kTH1F, {pAxis});
        histos.add("withtof/pt", "pt", kTH1F, {ptAxis});
        histos.add("withtof/length", "length", kTH1F, {lengthAxis});
        histos.add("withtof/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
        histos.add("withtof/beta", "beta", kTH2F, {pAxis, betaAxis});
        histos.add("withtof/delta", "delta", kTH2F, {pAxis, deltaAxis});
        histos.add("withtof/expP", "expP", kTH2F, {pAxis, pAxis});
        histos.add("withtof/mass", "mass", kTH1F, {massAxis});
        histos.add("withtof/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

        histos.addClone("withtof/", "goodreso/");
        histos.addClone("withtof/", "badreso/");
        histos.addClone("withtof/", "goodforevtime/");
        histos.addClone("withtof/", "withqualitycuts/");
      }

      listEfficiency.setObject(new THashList);
      auto makeEfficiency = [&](TString effname, TString efftitle) {
        listEfficiency->Add(new TEfficiency(effname, efftitle + ";TOF multiplicity;Efficiency", nBinsMultiplicity, 0, rangeMultiplicity));
      };

      makeEfficiency("effTOFEvTime", "Efficiency of the TOF Event Time");
      makeEfficiency("effT0ACEvTime", "Efficiency of the T0AC Event Time");
      makeEfficiency("effTOFT0ACEvTime", "Efficiency of the TOF+T0AC Event Time");
      makeEfficiency("effT0AEvTime", "Efficiency of the T0A Event Time");
      makeEfficiency("effT0CEvTime", "Efficiency of the T0C Event Time");
    }

    if (!doprocessMC) {
      return;
    }
    const AxisSpec diffAxis{1000, -1000, 1000, "Difference"};
    histos.add<TH1>("MC/diff/All", "All", HistType::kTH1F, {diffAxis})->GetXaxis()->SetTitle("t^{MC}_{ev}-t^{All}_{ev} (ps)");
    histos.add<TH1>("MC/diff/FT0", "FT0", HistType::kTH1F, {diffAxis})->GetXaxis()->SetTitle("t^{MC}_{ev}-t_{FT0}_{ev} (ps)");
    histos.add<TH1>("MC/diff/TOF", "TOF", HistType::kTH1F, {diffAxis})->GetXaxis()->SetTitle("t^{MC}_{ev}-t_{TOF}_{ev} (ps)");

    histos.add<TH2>("MC/diffvsZ/All", "All", HistType::kTH2F, {diffAxis, {100, -20, 20}})->GetXaxis()->SetTitle("t^{MC}_{ev}-t^{All}_{ev}_{ev} (ps)");
    histos.add<TH2>("MC/diffvsZ/FT0", "FT0", HistType::kTH2F, {diffAxis, {100, -20, 20}})->GetXaxis()->SetTitle("t^{MC}_{ev}-t_{FT0}_{ev} (ps)");
    histos.add<TH2>("MC/diffvsZ/TOF", "TOF", HistType::kTH2F, {diffAxis, {100, -20, 20}})->GetXaxis()->SetTitle("t^{MC}_{ev}-t_{TOF}_{ev} (ps)");

    // pion
    histos.add<TH2>("MC/particle/pdg211/all/particleDiff", "particleDiff", HistType::kTH2F, {ptAxis, diffAxis})->GetYaxis()->SetTitle("t^{MC}_{ev}-t^{part}_{MC} (ps)");
    histos.add<TH2>("MC/particle/pdg211/all/delta", "delta", HistType::kTH2F, {ptAxis, diffAxis})->GetYaxis()->SetTitle("t_{TOF}-t^{MC}_{ev}-t_{exp} (ps)");
    histos.add<TH2>("MC/particle/pdg211/all/deltaTRD", "deltaTRD", HistType::kTH2F, {ptAxis, diffAxis})->GetYaxis()->SetTitle("t_{TOF}-t^{MC}_{ev}-t_{exp} (ps)");
    histos.add<TH2>("MC/particle/pdg211/all/deltaNoTRD", "deltaNoTRD", HistType::kTH2F, {ptAxis, diffAxis})->GetYaxis()->SetTitle("t_{TOF}-t^{MC}_{ev}-t_{exp} (ps)");
    histos.addClone("MC/particle/pdg211/all/", "MC/particle/pdg211/prm/");

    histos.addClone("MC/particle/pdg211/", "MC/particle/pdgNeg211/");

    // kaon
    histos.addClone("MC/particle/pdg211/", "MC/particle/pdg321/");
    histos.addClone("MC/particle/pdg211/", "MC/particle/pdgNeg321/");
    // proton
    histos.addClone("MC/particle/pdg211/", "MC/particle/pdg2212/");
    histos.addClone("MC/particle/pdg211/", "MC/particle/pdgNeg2212/");

    const AxisSpec mcTimeAxis{1000, -1000, 1000, "MC coll time (ps)"};

    histos.add("MC/CollisionTime/eventtimeMC", "", HistType::kTH1F, {mcTimeAxis});
    histos.add<TH1>("MC/CollisionTime/All", "", HistType::kTH1F, {mcTimeAxis})->GetXaxis()->SetTitle("All (ps)");
    histos.add<TH1>("MC/CollisionTime/FT0", "", HistType::kTH1F, {mcTimeAxis})->GetXaxis()->SetTitle("FT0 (ps)");
    histos.add<TH1>("MC/CollisionTime/TOF", "", HistType::kTH1F, {mcTimeAxis})->GetXaxis()->SetTitle("TOF (ps)");
    histos.add<TH2>("MC/CollisionTime/eventtimeMCvsAll", "", HistType::kTH2F, {mcTimeAxis, mcTimeAxis})->GetYaxis()->SetTitle("All (ps)");
    histos.add<TH2>("MC/CollisionTime/eventtimeMCvsFT0", "", HistType::kTH2F, {mcTimeAxis, mcTimeAxis})->GetYaxis()->SetTitle("FT0 (ps)");
    histos.add<TH2>("MC/CollisionTime/eventtimeMCvsTOF", "", HistType::kTH2F, {mcTimeAxis, mcTimeAxis})->GetYaxis()->SetTitle("TOF (ps)");

    const AxisSpec axisBCID{o2::constants::lhc::LHCMaxBunches, -0.5, -0.5 + o2::constants::lhc::LHCMaxBunches, "BC ID in orbit"};
    const AxisSpec axisBCIDMC{o2::constants::lhc::LHCMaxBunches, -0.5, -0.5 + o2::constants::lhc::LHCMaxBunches, "MC BC ID in orbit"};

    histos.add("collisions/Reco/BCvsMCBC", "BC vs MC BC", kTH2D, {axisBCID, axisBCIDMC});
    histos.add<TH2>("collisions/Reco/FoundBCvsMCBC", "Found BC vs MC BC", kTH2D, {axisBCID, axisBCIDMC})->GetXaxis()->SetTitle("Found BC ID in orbit");
    histos.add<TH2>("collisions/Reco/FoundBCvsBC", "Found BC vs MC BC", kTH2D, {axisBCID, axisBCID})->GetXaxis()->SetTitle("Found BC ID in orbit");
    histos.add("collisions/Reco/bcMinusfoundBc", "bcMinusfoundBc", kTH1D, {{1600, -1000, 1000, "bc - foundBc (ns)"}});
    histos.add("collisions/Reco/bcMinusfoundBcRatio", "bcMinusfoundBcRatio", kTH1D, {{1600, -40, 40, "(bc - foundBc)/collisionTimeRes"}});
    histos.add("collisions/Reco/bcMinusMcBc", "bcMinusMcBc", kTH1D, {{1600, -1000, 1000, "bc - mcBc (ns)"}});
    histos.add("collisions/Reco/foundbcMinusMcBc", "foundbcMinusMcBc", kTH1D, {{1600, -1000, 1000, "foundBc - mcBc (ns)"}});
    histos.add("collisions/Reco/bcMinusMcBcRatio", "bcMinusMcBcRatio", kTH1D, {{1600, -40, 40, "(bc - mcBc)/collisionTimeRes"}});
    histos.add("collisions/Reco/foundbcMinusMcBcRatio", "foundbcMinusMcBcRatio", kTH1D, {{1600, -40, 40, "(foundBc-mcBc)/collisionTimeRes"}});

    histos.add("ft0/FT0AMinusFT0ACorrected", "", kTH1D, {{1600, -1000, 1000, "FT0A - FT0A_{corr} (ps)"}});
    histos.add("ft0/FT0CMinusFT0CCorrected", "", kTH1D, {{1600, -1000, 1000, "FT0C - FT0C_{corr} (ps)"}});
    histos.add("ft0/FT0ACMinusFT0ACCorrected", "", kTH1D, {{1600, -1000, 1000, "FT0AC - FT0AC_{corr} (ps)"}});
    histos.add("ft0/diffPosZ", "", kTH1D, {{100, -10, 10, "z_{mc} - z_{FT0} (cm)"}});
  }

  using Trks = soa::Join<aod::TracksIU, aod::TracksExtra,
                         aod::TOFSignal,
                         aod::TOFEvTime, aod::EvTimeTOFOnly, aod::pidEvTimeFlags,
                         aod::pidTOFbeta, aod::pidTOFmass>;
  using TrksData = soa::Join<Trks, aod::TrackSelection>;
  using EvTimeCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>;
  // Define slice per collision
  Preslice<Trks> perCollision = aod::track::collisionId;
  void processData(TrksData const& tracks, EvTimeCollisions const&)
  {
    static int ncolls = 0;
    int lastCollisionId = -1; // Last collision ID analysed
    for (auto& t : tracks) {
      if (!t.has_collision()) { // Track was not assigned to a collision
        continue;
      } else if (t.collisionId() == lastCollisionId) { // Event was already processed
        continue;
      }
      lastCollisionId = t.collisionId(); /// Cache last collision ID
      auto collision = t.collision_as<EvTimeCollisions>();

      histos.fill(HIST("eventSelection"), 0.5f);
      histos.fill(HIST("eventSelection"), 1.5f);
      if (t.tofEvTimeErr() > 199.f) {
        histos.fill(HIST("eventSelection"), 2.5f);
      } else {
        histos.fill(HIST("eventSelection"), 3.5f);
      }
      histos.fill(HIST("eventTime"), t.tofEvTime());
      histos.fill(HIST("eventTimeReso"), t.tofEvTimeErr());
      histos.fill(HIST("eventTimeVsMult"), t.evTimeTOFMult(), t.tofEvTime());
      histos.fill(HIST("eventTimeResoVsMult"), t.evTimeTOFMult(), t.tofEvTimeErr());

      if (t.isEvTimeTOF()) {
        histos.fill(HIST("eventTimeTOF"), t.evTimeTOF());
        histos.fill(HIST("eventTimeTOFReso"), t.evTimeTOFErr());
        histos.fill(HIST("eventTimeTOFMult"), t.evTimeTOFMult());
        histos.fill(HIST("eventTimeTOFVsMult"), t.evTimeTOFMult(), t.evTimeTOF());
        histos.fill(HIST("eventTimeTOFResoVsMult"), t.evTimeTOFMult(), t.evTimeTOFErr());
      }

      if (collision.has_foundFT0()) { // T0 measurement is available
        if (collision.t0ACorrectedValid()) {
          histos.fill(HIST("eventTimeT0A"), collision.t0ACorrected() * 1000.f);
          if (t.isEvTimeTOF()) {
            histos.fill(HIST("deltaEvTimeTOFT0A"), t.evTimeTOF() - collision.t0ACorrected() * 1000.f);
            histos.fill(HIST("deltaEvTimeTOFT0AvsTOF"), t.evTimeTOF(), t.evTimeTOF() - collision.t0ACorrected() * 1000.f);
          }
        }
        if (collision.t0CCorrectedValid()) {
          histos.fill(HIST("eventTimeT0C"), collision.t0CCorrected() * 1000.f);
          if (t.isEvTimeTOF()) {
            histos.fill(HIST("deltaEvTimeTOFT0C"), t.evTimeTOF() - collision.t0CCorrected() * 1000.f);
            histos.fill(HIST("deltaEvTimeTOFT0CvsTOF"), t.evTimeTOF(), t.evTimeTOF() - collision.t0CCorrected() * 1000.f);
          }
        }
        if (collision.t0ACValid()) {
          histos.fill(HIST("eventTimeT0AC"), collision.t0AC() * 1000.f);
          histos.fill(HIST("eventTimeT0ACReso"), collision.t0resolution() * 1000.f);
          if (t.isEvTimeTOF()) {
            histos.fill(HIST("deltaEvTimeTOFT0AC"), t.evTimeTOF() - collision.t0AC() * 1000.f);
            histos.fill(HIST("deltaEvTimeTOFT0ACvsTOF"), t.evTimeTOF(), t.evTimeTOF() - collision.t0AC() * 1000.f);
          }
        }
        if (collision.t0ACorrectedValid() && collision.t0CCorrectedValid() && t.isEvTimeTOF()) {
          histos.fill(HIST("deltaEvTimeTOFT0AvsT0C"), t.evTimeTOF() - collision.t0ACorrected() * 1000.f, t.evTimeTOF() - collision.t0CCorrected() * 1000.f);
        }
      }

      histos.fill(HIST("collisionTime"), collision.collisionTime());
      histos.fill(HIST("collisionTimeRes"), collision.collisionTimeRes());
      ncolls++;

      const auto tracksInCollision = tracks.sliceBy(perCollision, lastCollisionId);
      int nTracksWithTOF = 0;
      for (auto const& trk : tracksInCollision) { // Loop on Tracks
        histos.fill(HIST("trackSelection"), 0.5f);

        if (!trk.isGlobalTrack()) {
          continue;
        }
        histos.fill(HIST("trackSelection"), 1.5f);

        if (!trk.hasITS()) {
          continue;
        }
        histos.fill(HIST("trackSelection"), 2.5f);
        if (!trk.hasTPC()) {
          continue;
        }
        histos.fill(HIST("trackSelection"), 3.5f);

        histos.fill(HIST("tracks/p"), trk.p());
        histos.fill(HIST("tracks/pt"), trk.pt());
        histos.fill(HIST("tracks/length"), trk.length());

        if (enableDebug) {
          if (trk.tofEvTimeErr() > 199.f) {
            histos.fill(HIST("badreso/ptden"), trk.pt());
          } else {
            histos.fill(HIST("goodreso/ptden"), trk.pt());
          }
        }

        if (!trk.hasTOF()) {
          continue;
        }
        nTracksWithTOF++;
        histos.fill(HIST("trackSelection"), 4.5f);

        // Recompute quantities with event times
        const float& betaTOF = trk.evTimeTOFMult() > 1 ? o2::pid::tof::Beta::GetBeta(trk, trk.evTimeTOF()) : 999.f;
        const float& betaT0A = collision.t0ACorrectedValid() ? o2::pid::tof::Beta::GetBeta(trk, collision.t0ACorrected() * 1000.f) : 999.f;
        const float& betaT0C = collision.t0CCorrectedValid() ? o2::pid::tof::Beta::GetBeta(trk, collision.t0CCorrected() * 1000.f) : 999.f;
        const float& betaT0AC = collision.t0ACValid() ? o2::pid::tof::Beta::GetBeta(trk, collision.t0AC() * 1000.f) : 999.f;

        const float& massTOF = trk.evTimeTOFMult() > 1 ? o2::pid::tof::TOFMass::GetTOFMass(trk, betaTOF) : 999.f;
        const float& massT0A = collision.t0ACorrectedValid() ? o2::pid::tof::TOFMass::GetTOFMass(trk, betaT0A) : 999.f;
        const float& massT0C = collision.t0CCorrectedValid() ? o2::pid::tof::TOFMass::GetTOFMass(trk, betaT0C) : 999.f;
        const float& massT0AC = collision.t0ACValid() ? o2::pid::tof::TOFMass::GetTOFMass(trk, betaT0AC) : 999.f;

        const float& deltaPi = trk.tofSignal() - trk.tofEvTime() - trk.tofExpTimePi();

        histos.fill(HIST("tofbeta/inclusive"), trk.p(), trk.beta());
        histos.fill(HIST("tofmass/inclusive"), trk.p(), trk.mass());
        if (trk.isEvTimeTOF()) {
          histos.fill(HIST("tofbeta/EvTimeTOF"), trk.p(), trk.beta());
          histos.fill(HIST("tofmass/EvTimeTOF"), trk.p(), trk.mass());
          histos.fill(HIST("tofbeta/EvTimeTOFOnly"), trk.p(), betaTOF);
          histos.fill(HIST("tofmass/EvTimeTOFOnly"), trk.p(), massTOF);
        }
        if (trk.isEvTimeT0AC()) {
          histos.fill(HIST("tofbeta/EvTimeT0AC"), trk.p(), trk.beta());
          histos.fill(HIST("tofmass/EvTimeT0AC"), trk.p(), trk.mass());
        }
        histos.fill(HIST("tofbeta/EvTimeT0AOnly"), trk.p(), betaT0A);
        histos.fill(HIST("tofbeta/EvTimeT0COnly"), trk.p(), betaT0C);
        histos.fill(HIST("tofbeta/EvTimeT0ACOnly"), trk.p(), betaT0AC);

        histos.fill(HIST("tofmass/EvTimeT0AOnly"), trk.p(), massT0A);
        histos.fill(HIST("tofmass/EvTimeT0COnly"), trk.p(), massT0C);
        histos.fill(HIST("tofmass/EvTimeT0ACOnly"), trk.p(), massT0AC);

        if (trk.p() > minPReso && trk.p() < maxPReso) {
          histos.fill(HIST("deltaVsMult/pi"), trk.evTimeTOFMult(), deltaPi);
          histos.fill(HIST("deltaVsReso/pi"), trk.evTimeTOFMult(), deltaPi);
        }
        if (enableDebug) {

          histos.fill(HIST("withtof/p"), trk.p());
          histos.fill(HIST("withtof/pt"), trk.pt());
          histos.fill(HIST("withtof/length"), trk.length());
          histos.fill(HIST("withtof/tofSignal"), trk.tofSignal());
          histos.fill(HIST("withtof/beta"), trk.p(), trk.beta());
          histos.fill(HIST("withtof/delta"), trk.p(), deltaPi);

          histos.fill(HIST("withtof/expP"), trk.p(), trk.tofExpMom());
          histos.fill(HIST("withtof/mass"), trk.mass());
          histos.fill(HIST("withtof/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
          if (trk.pt() > 0.3 && trk.beta() > 0.3) {
            histos.fill(HIST("withqualitycuts/p"), trk.p());
            histos.fill(HIST("withqualitycuts/pt"), trk.pt());
            histos.fill(HIST("withqualitycuts/length"), trk.length());
            histos.fill(HIST("withqualitycuts/mass"), trk.mass());
          }

          if (trk.tofEvTimeErr() > 199.f) {
            histos.fill(HIST("badreso/p"), trk.p());
            histos.fill(HIST("badreso/pt"), trk.pt());
            histos.fill(HIST("badreso/length"), trk.length());
            histos.fill(HIST("badreso/tofSignal"), trk.tofSignal());
            histos.fill(HIST("badreso/beta"), trk.p(), trk.beta());
            histos.fill(HIST("badreso/delta"), trk.p(), deltaPi);
            histos.fill(HIST("badreso/expP"), trk.p(), trk.tofExpMom());
            histos.fill(HIST("badreso/mass"), trk.mass());
            histos.fill(HIST("badreso/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
          } else {
            histos.fill(HIST("goodreso/p"), trk.p());
            histos.fill(HIST("goodreso/pt"), trk.pt());
            histos.fill(HIST("goodreso/length"), trk.length());
            histos.fill(HIST("goodreso/tofSignal"), trk.tofSignal());
            histos.fill(HIST("goodreso/beta"), trk.p(), trk.beta());
            histos.fill(HIST("goodreso/delta"), trk.p(), deltaPi);
            histos.fill(HIST("goodreso/expP"), trk.p(), trk.tofExpMom());
            histos.fill(HIST("goodreso/mass"), trk.mass());
            histos.fill(HIST("goodreso/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
          }
          if (!trk.usedForTOFEvTime()) {
            continue;
          }
          histos.fill(HIST("goodforevtime/p"), trk.p());
          histos.fill(HIST("goodforevtime/pt"), trk.pt());
          histos.fill(HIST("goodforevtime/length"), trk.length());
          histos.fill(HIST("goodforevtime/tofSignal"), trk.tofSignal());
          histos.fill(HIST("goodforevtime/beta"), trk.p(), trk.beta());
          histos.fill(HIST("goodforevtime/delta"), trk.p(), deltaPi);
          histos.fill(HIST("goodforevtime/expP"), trk.p(), trk.tofExpMom());
          histos.fill(HIST("goodforevtime/mass"), trk.mass());
          histos.fill(HIST("goodforevtime/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
        }
      }
      static_cast<TEfficiency*>(listEfficiency->FindObject("effTOFEvTime"))->Fill(t.isEvTimeTOF(), nTracksWithTOF);
      static_cast<TEfficiency*>(listEfficiency->FindObject("effT0AEvTime"))->Fill(collision.has_foundFT0() && collision.t0ACorrectedValid(), nTracksWithTOF);
      static_cast<TEfficiency*>(listEfficiency->FindObject("effT0CEvTime"))->Fill(collision.has_foundFT0() && collision.t0CCorrectedValid(), nTracksWithTOF);
      static_cast<TEfficiency*>(listEfficiency->FindObject("effT0ACEvTime"))->Fill(collision.has_foundFT0() && collision.t0ACValid(), nTracksWithTOF);
      static_cast<TEfficiency*>(listEfficiency->FindObject("effTOFT0ACEvTime"))->Fill(t.isEvTimeTOF() && collision.has_foundFT0() && collision.t0ACorrectedValid(), nTracksWithTOF);
    }
  }
  PROCESS_SWITCH(tofPidCollisionTimeQa, processData, "Process data", true);

  using TrksMC = soa::Join<Trks, aod::McTrackLabels>;
  using EvTimeCollisionsMC = soa::Join<EvTimeCollisions, aod::McCollisionLabels>;
  void processMC(TrksMC const& tracks,
                 EvTimeCollisionsMC const&,
                 aod::McParticles const&,
                 aod::BCs const&,
                 aod::FT0s const&,
                 aod::McCollisions const&)
  {
    // static int ncolls = 0;
    int lastCollisionId = -1; // Last collision ID analysed

    for (auto& trk : tracks) {
      if (!trk.has_collision()) { // Track was not assigned to a collision
        continue;
      }
      const auto& collision = trk.collision_as<EvTimeCollisionsMC>();
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& collisionMC = collision.mcCollision_as<aod::McCollisions>();
      const float eventtimeMC = collisionMC.t() * 1000.f;

      if (trk.has_mcParticle()) {
        const auto& particle = trk.mcParticle();
        const auto& mcCollTimeMinusFormationTime = particle.vt() - collisionMC.t();
        const auto& mcTOFvalue = trk.tofSignal() - eventtimeMC - trk.tofExpTimePi();
        LOG(debug) << "Track " << particle.vt() << " vs " << eventtimeMC;
        switch (particle.pdgCode()) {
          case 211:
            histos.fill(HIST("MC/particle/pdg211/all/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            histos.fill(HIST("MC/particle/pdg211/all/delta"), particle.pt(), mcTOFvalue);
            if (trk.hasTRD()) {
              histos.fill(HIST("MC/particle/pdg211/all/deltaTRD"), particle.pt(), mcTOFvalue);
            } else {
              histos.fill(HIST("MC/particle/pdg211/all/deltaNoTRD"), particle.pt(), mcTOFvalue);
            }
            if (particle.isPhysicalPrimary()) {
              histos.fill(HIST("MC/particle/pdg211/prm/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            }
            break;
          case -211:
            histos.fill(HIST("MC/particle/pdgNeg211/all/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            histos.fill(HIST("MC/particle/pdgNeg211/all/delta"), particle.pt(), mcTOFvalue);
            histos.fill(HIST("MC/particle/pdgNeg211/all/delta"), particle.pt(), mcTOFvalue);
            if (trk.hasTRD()) {
              histos.fill(HIST("MC/particle/pdgNeg211/all/deltaTRD"), particle.pt(), mcTOFvalue);
            } else {
              histos.fill(HIST("MC/particle/pdgNeg211/all/deltaNoTRD"), particle.pt(), mcTOFvalue);
            }
            if (particle.isPhysicalPrimary()) {
              histos.fill(HIST("MC/particle/pdgNeg211/prm/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            }
            break;
          case 321:
            histos.fill(HIST("MC/particle/pdg321/all/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            if (particle.isPhysicalPrimary()) {
              histos.fill(HIST("MC/particle/pdg321/prm/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            }
            break;
          case -321:
            histos.fill(HIST("MC/particle/pdgNeg321/all/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            if (particle.isPhysicalPrimary()) {
              histos.fill(HIST("MC/particle/pdgNeg321/prm/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            }
            break;
          case 2212:
            histos.fill(HIST("MC/particle/pdg2212/all/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            if (particle.isPhysicalPrimary()) {
              histos.fill(HIST("MC/particle/pdg2212/prm/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            }
            break;
          case -2212:
            histos.fill(HIST("MC/particle/pdgNeg2212/all/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            if (particle.isPhysicalPrimary()) {
              histos.fill(HIST("MC/particle/pdgNeg2212/prm/particleDiff"), particle.pt(), mcCollTimeMinusFormationTime);
            }
            break;
          default:
            break;
        }
      }

      if (trk.collisionId() == lastCollisionId) { // Event was already processed
        continue;
      }
      lastCollisionId = trk.collisionId(); /// Cache last collision ID

      float t0AC[2] = {0.f, 0.f};
      // const auto& ft0 = collision.foundFT0();
      if (collision.t0ACValid()) {
        t0AC[0] = collision.t0AC() * 1000.f;
        t0AC[1] = collision.t0resolution() * 1000.f;
      }

      float t0A = 1e10;
      float t0C = 1e10;

      constexpr float dummyTime = 30.; // Due to HW limitations time can be only within range (-25,25) ns, dummy time is around 32 ns
      if (collision.has_foundFT0()) {  // Get the non corrected FT0AC
        const auto& ft0 = collision.foundFT0();
        const std::bitset<8>& triggers = ft0.triggerMask();
        const bool ora = triggers[o2::ft0::Triggers::bitA];
        const bool orc = triggers[o2::ft0::Triggers::bitC];
        if (ora && ft0.timeA() < dummyTime) {
          t0A = ft0.timeA();
        }
        if (orc && ft0.timeC() < dummyTime) {
          t0C = ft0.timeC();
        }
        if (ft0.isValidTime()) {
          histos.fill(HIST("ft0/diffPosZ"), ft0.posZ() - collisionMC.posZ());
        }
      }
      if (collision.t0ACorrectedValid()) {
        histos.fill(HIST("ft0/FT0AMinusFT0ACorrected"), (t0A - collision.t0ACorrected()) * 1000.f);
      }
      if (collision.t0CCorrectedValid()) {
        histos.fill(HIST("ft0/FT0CMinusFT0CCorrected"), (t0C - collision.t0CCorrected()) * 1000.f);
      }
      if (collision.t0ACValid()) {
        histos.fill(HIST("ft0/FT0ACMinusFT0ACCorrected"), (0.5 * (t0A + t0C) - collision.t0AC()) * 1000.f);
      }

      const auto& recoBC = collision.bc();
      const auto& foundBC = collision.foundBC();
      const auto& mcBC = collisionMC.bc();

      const auto& recoBCid = recoBC.globalBC() % o2::constants::lhc::LHCMaxBunches;
      const auto& mcBCid = mcBC.globalBC() % o2::constants::lhc::LHCMaxBunches;
      const auto& foundBCid = foundBC.globalBC() % o2::constants::lhc::LHCMaxBunches;
      const int diffRecoFoundBC = foundBC.globalBC() - recoBC.globalBC();
      const int diffRecoMCBC = recoBC.globalBC() - mcBC.globalBC();
      const int diffFoundMCBC = foundBC.globalBC() - mcBC.globalBC();
      histos.fill(HIST("collisions/Reco/BCvsMCBC"), recoBCid, mcBCid);
      histos.fill(HIST("collisions/Reco/FoundBCvsMCBC"), foundBCid, mcBCid);
      histos.fill(HIST("collisions/Reco/FoundBCvsBC"), foundBCid, recoBCid);
      histos.fill(HIST("collisions/Reco/bcMinusMcBc"), (diffRecoMCBC)*o2::constants::lhc::LHCBunchSpacingNS);
      histos.fill(HIST("collisions/Reco/foundbcMinusMcBc"), (diffFoundMCBC)*o2::constants::lhc::LHCBunchSpacingNS);
      histos.fill(HIST("collisions/Reco/bcMinusfoundBc"), (diffRecoFoundBC)*o2::constants::lhc::LHCBunchSpacingNS);
      histos.fill(HIST("collisions/Reco/bcMinusfoundBcRatio"), (diffRecoFoundBC)*o2::constants::lhc::LHCBunchSpacingNS / collision.collisionTimeRes());
      histos.fill(HIST("collisions/Reco/foundbcMinusMcBcRatio"), (diffFoundMCBC)*o2::constants::lhc::LHCBunchSpacingNS / collision.collisionTimeRes());
      histos.fill(HIST("collisions/Reco/bcMinusMcBcRatio"), (diffRecoMCBC)*o2::constants::lhc::LHCBunchSpacingNS / collision.collisionTimeRes());

      // timeInBCNS + bc2ns();
      // bc* o2::constants::lhc::LHCBunchSpacingNS + orbit* o2::constants::lhc::LHCOrbitNS;
      // int64_t(mcBC.globalBC() * o2::constants::lhc::LHCBunchSpacingNS * 1e-3))

      histos.fill(HIST("MC/CollisionTime/eventtimeMC"), eventtimeMC);
      histos.fill(HIST("MC/CollisionTime/All"), trk.tofEvTime());
      histos.fill(HIST("MC/CollisionTime/FT0"), t0AC[0]);
      histos.fill(HIST("MC/CollisionTime/TOF"), trk.evTimeTOF());
      histos.fill(HIST("MC/CollisionTime/eventtimeMCvsAll"), eventtimeMC, trk.tofEvTime());
      histos.fill(HIST("MC/CollisionTime/eventtimeMCvsFT0"), eventtimeMC, t0AC[0]);
      histos.fill(HIST("MC/CollisionTime/eventtimeMCvsTOF"), eventtimeMC, trk.evTimeTOF());

      histos.fill(HIST("MC/diff/All"), eventtimeMC - trk.tofEvTime());
      histos.fill(HIST("MC/diff/FT0"), eventtimeMC - t0AC[0]);
      histos.fill(HIST("MC/diff/TOF"), eventtimeMC - trk.evTimeTOF());
      histos.fill(HIST("MC/diffvsZ/All"), eventtimeMC - trk.tofEvTime(), collisionMC.posZ());
      histos.fill(HIST("MC/diffvsZ/FT0"), eventtimeMC - t0AC[0], collisionMC.posZ());
      histos.fill(HIST("MC/diffvsZ/TOF"), eventtimeMC - trk.evTimeTOF(), collisionMC.posZ());
    }
  }
  PROCESS_SWITCH(tofPidCollisionTimeQa, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tofPidCollisionTimeQa>(cfgc)}; }
