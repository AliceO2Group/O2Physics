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
/// \file   pidTOFFull.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce PID tables for TOF split for each particle.
///         The event time maker can be used to produce event TOF times.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "TableHelper.h"
#include "Framework/StaticFor.h"
#include "TOFBase/EventTimeMaker.h"
#include "pidTOFBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TOF PID QA histograms"}}};
  options.push_back({"add-qa-ev-time", VariantType::Int, 0, {"Produce TOF PID QA histograms for TOF event time"}});
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Table with the TOF event time
namespace o2::aod
{
namespace tofeventtime
{
DECLARE_SOA_COLUMN(TOFEvTime, tofEvTime, float);       //! TOF event time
DECLARE_SOA_COLUMN(TOFEvTimeErr, tofEvTimeErr, float); //! TOF event time error
DECLARE_SOA_COLUMN(TOFEvTimeMult, tofEvTimeMult, int); //! TOF event time multiplicity
} // namespace tofeventtime

DECLARE_SOA_TABLE(TOFEvTime, "AOD", "TOFEvTime", //! Table of the TOF event time
                  tofeventtime::TOFEvTime,
                  tofeventtime::TOFEvTimeErr,
                  tofeventtime::TOFEvTimeMult);
} // namespace o2::aod

/// Task to produce the response table
struct tofPidFull {
  // Tables to produce
  Produces<o2::aod::TOFEvTime> tableEvTime;
  Produces<o2::aod::pidTOFFullEl> tablePIDEl;
  Produces<o2::aod::pidTOFFullMu> tablePIDMu;
  Produces<o2::aod::pidTOFFullPi> tablePIDPi;
  Produces<o2::aod::pidTOFFullKa> tablePIDKa;
  Produces<o2::aod::pidTOFFullPr> tablePIDPr;
  Produces<o2::aod::pidTOFFullDe> tablePIDDe;
  Produces<o2::aod::pidTOFFullTr> tablePIDTr;
  Produces<o2::aod::pidTOFFullHe> tablePIDHe;
  Produces<o2::aod::pidTOFFullAl> tablePIDAl;
  // Detector response and input parameters
  DetectorResponse response;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> sigmaname{"param-sigma", "TOFReso", "Name of the parametrization for the expected sigma, used in both file and CCDB mode"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TOF", "Path of the TOF parametrization on the CCDB"};
  Configurable<long> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  // Configuration flags to include and exclude particle hypotheses
  Configurable<int> pidEl{"pid-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidMu{"pid-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPi{"pid-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidKa{"pid-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPr{"pid-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidDe{"pid-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTr{"pid-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidHe{"pid-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidAl{"pid-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};

  void init(o2::framework::InitContext& initContext)
  {
    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      const std::string table = "pidTOFFull" + particle;
      if (isTableRequiredInWorkflow(initContext, table)) {
        if (flag < 0) {
          flag.value = 1;
          LOG(info) << "Auto-enabling table: " + table;
        } else if (flag > 0) {
          flag.value = 1;
          LOG(info) << "Table enabled: " + table;
        } else {
          LOG(info) << "Table disabled: " + table;
        }
      }
    };

    enableFlag("El", pidEl);
    enableFlag("Mu", pidMu);
    enableFlag("Pi", pidPi);
    enableFlag("Ka", pidKa);
    enableFlag("Pr", pidPr);
    enableFlag("De", pidDe);
    enableFlag("Tr", pidTr);
    enableFlag("He", pidHe);
    enableFlag("Al", pidAl);

    // Getting the parametrization parameters
    ccdb->setURL(url.value);
    ccdb->setTimestamp(timestamp.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    //
    const std::vector<float> p = {0.008, 0.008, 0.002, 40.0};
    response.SetParameters(DetectorResponse::kSigma, p);
    const std::string fname = paramfile.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file" << fname << ", using param: " << sigmaname.value;
      response.LoadParamFromFile(fname.data(), sigmaname.value, DetectorResponse::kSigma);
    } else { // Loading it from CCDB
      std::string path = ccdbPath.value + "/" + sigmaname.value;
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << path << " for timestamp " << timestamp.value;
      response.LoadParam(DetectorResponse::kSigma, ccdb->getForTimeStamp<Parametrization>(path, timestamp.value));
    }
  }

  using TrksEvTime = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection>;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksEvTime::iterator, pid>;
  void processEvTime(TrksEvTime const& tracks, aod::Collisions const&)
  {
    constexpr auto responseEl = ResponseImplementationEvTime<PID::Electron>();
    constexpr auto responseMu = ResponseImplementationEvTime<PID::Muon>();
    constexpr auto responsePi = ResponseImplementationEvTime<PID::Pion>();
    constexpr auto responseKa = ResponseImplementationEvTime<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementationEvTime<PID::Proton>();
    constexpr auto responseDe = ResponseImplementationEvTime<PID::Deuteron>();
    constexpr auto responseTr = ResponseImplementationEvTime<PID::Triton>();
    constexpr auto responseHe = ResponseImplementationEvTime<PID::Helium3>();
    constexpr auto responseAl = ResponseImplementationEvTime<PID::Alpha>();

    tableEvTime.reserve(tracks.size());

    auto reserveTable = [&tracks](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(tracks.size());
    };

    reserveTable(pidEl, tablePIDEl);
    reserveTable(pidMu, tablePIDMu);
    reserveTable(pidPi, tablePIDPi);
    reserveTable(pidKa, tablePIDKa);
    reserveTable(pidPr, tablePIDPr);
    reserveTable(pidDe, tablePIDDe);
    reserveTable(pidTr, tablePIDTr);
    reserveTable(pidHe, tablePIDHe);
    reserveTable(pidAl, tablePIDAl);

    int lastCollisionId = -1;      // Last collision ID analysed
    for (auto const& t : tracks) { // Loop on collisions
      if (!t.has_collision()) {    // Track was not assigned, cannot compute event time
        tableEvTime(0.f, 999.f, -1);

        auto fillEmptyTable = [](const Configurable<int>& flag, auto& table) {
          if (flag.value != 1) {
            return;
          }
          table(-999.f, -999.f);
        };

        fillEmptyTable(pidEl, tablePIDEl);
        fillEmptyTable(pidMu, tablePIDMu);
        fillEmptyTable(pidPi, tablePIDPi);
        fillEmptyTable(pidKa, tablePIDKa);
        fillEmptyTable(pidPr, tablePIDPr);
        fillEmptyTable(pidDe, tablePIDDe);
        fillEmptyTable(pidTr, tablePIDTr);
        fillEmptyTable(pidHe, tablePIDHe);
        fillEmptyTable(pidAl, tablePIDAl);

        continue;
      } else if (t.collisionId() == lastCollisionId) { // Event time from this collision is already in the table
        continue;
      }
      /// Create new table for the tracks in a collision
      lastCollisionId = t.collisionId(); /// Cache last collision ID

      const auto tracksInCollision = tracks.sliceBy(aod::track::collisionId, lastCollisionId);
      // First make table for event time
      const auto evTime = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, response);
      static constexpr bool removebias = true;
      int ngoodtracks = 0;
      for (auto const& trk : tracksInCollision) { // Loop on Tracks
        float et = evTime.mEventTime;
        float erret = evTime.mEventTimeError;
        if constexpr (removebias) {
          evTime.removeBias<TrksEvTime::iterator, filterForTOFEventTime>(trk, ngoodtracks, et, erret);
        }
        tableEvTime(et, erret, evTime.mEventTimeMultiplicity);
      }

      // Check and fill enabled tables
      auto makeTable = [&tracksInCollision, &evTime, &ngoodtracks, this](const Configurable<int>& flag, auto& table, const auto& responsePID) {
        if (flag.value == 1) {
          ngoodtracks = 0;
          // Prepare memory for enabled tables
          table.reserve(tracksInCollision.size());
          for (auto const& trk : tracksInCollision) { // Loop on Tracks

            float et = evTime.mEventTime;
            float erret = evTime.mEventTimeError;
            if constexpr (removebias) {
              evTime.removeBias<TrksEvTime::iterator, filterForTOFEventTime>(trk, ngoodtracks, et, erret);
            }

            if (erret > 199.f) {
              table(erret, 999.f);
              continue;
            }

            table(responsePID.GetExpectedSigma(response, trk, trk.tofSignal(), erret),
                  responsePID.GetSeparation(response, trk, et, erret));
          }
        }
      };

      makeTable(pidEl, tablePIDEl, responseEl);
      makeTable(pidMu, tablePIDMu, responseMu);
      makeTable(pidPi, tablePIDPi, responsePi);
      makeTable(pidKa, tablePIDKa, responseKa);
      makeTable(pidPr, tablePIDPr, responsePr);
      makeTable(pidDe, tablePIDDe, responseDe);
      makeTable(pidTr, tablePIDTr, responseTr);
      makeTable(pidHe, tablePIDHe, responseHe);
      makeTable(pidAl, tablePIDAl, responseAl);
    }
  }

  PROCESS_SWITCH(tofPidFull, processEvTime, "Produce TOF response with TOF event time", false);

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal>;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Trks::iterator, pid>;
  void processNoEvTime(Trks const& tracks, aod::Collisions const&)
  {
    constexpr auto responseEl = ResponseImplementation<PID::Electron>();
    constexpr auto responseMu = ResponseImplementation<PID::Muon>();
    constexpr auto responsePi = ResponseImplementation<PID::Pion>();
    constexpr auto responseKa = ResponseImplementation<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementation<PID::Proton>();
    constexpr auto responseDe = ResponseImplementation<PID::Deuteron>();
    constexpr auto responseTr = ResponseImplementation<PID::Triton>();
    constexpr auto responseHe = ResponseImplementation<PID::Helium3>();
    constexpr auto responseAl = ResponseImplementation<PID::Alpha>();

    // Check and fill enabled tables
    auto makeTable = [&tracks](const Configurable<int>& flag, auto& table, const DetectorResponse& response, const auto& responsePID) {
      if (flag.value == 1) {
        // Prepare memory for enabled tables
        table.reserve(tracks.size());
        for (auto const& trk : tracks) { // Loop on Tracks
          table(responsePID.GetExpectedSigma(response, trk),
                responsePID.GetSeparation(response, trk));
        }
      }
    };
    makeTable(pidEl, tablePIDEl, response, responseEl);
    makeTable(pidMu, tablePIDMu, response, responseMu);
    makeTable(pidPi, tablePIDPi, response, responsePi);
    makeTable(pidKa, tablePIDKa, response, responseKa);
    makeTable(pidPr, tablePIDPr, response, responsePr);
    makeTable(pidDe, tablePIDDe, response, responseDe);
    makeTable(pidTr, tablePIDTr, response, responseTr);
    makeTable(pidHe, tablePIDHe, response, responseHe);
    makeTable(pidAl, tablePIDAl, response, responseAl);
  }

  PROCESS_SWITCH(tofPidFull, processNoEvTime, "Produce TOF response without TOF event time, standard for Run 2", true);
};

struct tofPidCollisionTimeQa { /// Task that checks the TOF collision time
  Configurable<int> nBinsEvTime{"nBinsEvTime", 1000, "Number of bins for the event time"};
  Configurable<float> minEvTime{"minEvTime", -1000.f, "Minimum in range in event time"};
  Configurable<float> maxEvTime{"maxEvTime", 1000.f, "Maximum in range in event time"};
  Configurable<int> nBinsTofSignal{"nBinsTofSignal", 5000, "Number of bins for the tof signal time"};
  Configurable<float> minTofSignal{"minTofSignal", 0.f, "Minimum in range in tof signal time"};
  Configurable<float> maxTofSignal{"maxTofSignal", 100e3, "Maximum in range in tof signal time"};
  Configurable<float> rangeEvTimeReso{"rangeEvTimeReso", 1000.f, "Range in event time resolution"};
  Configurable<int> nBinsMultiplicity{"nBinsMultiplicity", 1000, "Number of bins for the multiplicity"};
  Configurable<float> rangeMultiplicity{"rangeMultiplicity", 1000.f, "Range in event time resolution"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};
  void init(o2::framework::InitContext& initContext)
  {
    const AxisSpec evTimeAxis{nBinsEvTime, minEvTime, maxEvTime, "TOF event time (ps)"};
    const AxisSpec multAxis{nBinsEvTime, 0, rangeMultiplicity, "Track multiplicity for TOF event time"};
    const AxisSpec evTimeResoAxis{nBinsMultiplicity, 0, rangeEvTimeReso, "TOF event time resolution (ps)"};
    const AxisSpec tofSignalAxis{nBinsTofSignal, minTofSignal, maxTofSignal, "TOF signal (ps)"};
    AxisSpec pAxis{1000, 0.01, 5, "#it{p} GeV/#it{c}"};
    pAxis.makeLogaritmic();
    const AxisSpec ptAxis{100, 0, 10, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec collisionAxis{6000, -0.5f, 6000.f - .5f, "Collision index % 6000"};
    const AxisSpec massAxis{1000, 0, 3, "TOF mass (GeV/#it{c}^{2})"};
    const AxisSpec betaAxis{1000, 0, 1.5, "TOF #beta"};
    const AxisSpec deltaAxis{1000, -10000, 10000, "t-texp-t0"};
    const AxisSpec lengthAxis{1000, 0, 600, "Track length (cm)"};

    histos.add("eventSelection", "eventSelection", kTH1F, {{10, 0, 10}});
    histos.add("eventTime", "eventTime", kTH1F, {evTimeAxis});
    histos.add("eventTimeM", "eventTimeM", kTH1F, {evTimeAxis});
    histos.add("eventTimeReso", "eventTimeReso", kTH1F, {evTimeResoAxis});
    histos.add("eventTimeMult", "eventTimeMult", kTH1F, {multAxis});
    histos.add("collisionTime", "collisionTime", kTH1F, {evTimeResoAxis});
    histos.add("collisionTimeRes", "collisionTimeRes", kTH1F, {evTimeResoAxis});

    histos.add("tracks/p", "p", kTH1F, {pAxis});
    histos.add("tracks/pt", "pt", kTH1F, {ptAxis});
    histos.add("tracks/length", "length", kTH1F, {lengthAxis});

    histos.add("withtof/p", "p", kTH1F, {pAxis});
    histos.add("withtof/pt", "pt", kTH1F, {ptAxis});
    histos.add("withtof/length", "length", kTH1F, {lengthAxis});
    histos.add("withtof/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("withtof/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("withtof/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("withtof/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("withtof/mass", "mass", kTH1F, {massAxis});
    histos.add("withtof/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("goodforevtime/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("goodforevtime/p", "p", kTH1F, {pAxis});
    histos.add("goodforevtime/pt", "pt", kTH1F, {ptAxis});
    histos.add("goodforevtime/length", "length", kTH1F, {lengthAxis});
    histos.add("goodforevtime/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("goodforevtime/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("goodforevtime/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("goodforevtime/mass", "mass", kTH1F, {massAxis});
    histos.add("goodforevtime/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection>;
  int ncolls = 0;
  void process(aod::Collision const&, Trks const& tracks)
  {
    histos.fill(HIST("eventSelection"), 0.5f);
    bool eventSet = false;
    for (auto& t : tracks) {
      if (!t.isGlobalTrack()) {
        continue;
      }
      histos.fill(HIST("tracks/p"), t.p());
      histos.fill(HIST("tracks/pt"), t.pt());
      histos.fill(HIST("tracks/length"), t.length());
      if (!t.hasTOF()) {
        continue;
      }
      const float beta = o2::pid::tof::Beta<Trks::iterator>::GetBeta(t, t.tofEvTime());
      const float mass = o2::pid::tof::TOFMass<Trks::iterator>::GetTOFMass(t.p(), beta);
      histos.fill(HIST("withtof/p"), t.p());
      histos.fill(HIST("withtof/pt"), t.pt());
      histos.fill(HIST("withtof/length"), t.length());
      histos.fill(HIST("withtof/tofSignal"), t.tofSignal());
      histos.fill(HIST("withtof/beta"), t.p(), beta);
      histos.fill(HIST("withtof/delta"), t.p(), t.tofSignal() - t.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(t));
      histos.fill(HIST("withtof/expP"), t.p(), t.tofExpMom());
      histos.fill(HIST("withtof/mass"), mass);
      histos.fill(HIST("withtof/tofSignalPerCollision"), ncolls % 6000, t.tofSignal());
      if (!eventSet) {
        histos.fill(HIST("eventTime"), t.tofEvTime());
        if (t.tofEvTimeMult() > 1) {
          histos.fill(HIST("eventTimeM"), t.tofEvTime());
        }
        histos.fill(HIST("eventTimeReso"), t.tofEvTimeErr());
        histos.fill(HIST("eventTimeMult"), t.tofEvTimeMult());
        histos.fill(HIST("collisionTime"), t.collision().collisionTime());
        histos.fill(HIST("collisionTimeRes"), t.collision().collisionTimeRes());
        eventSet = true;
        ncolls++;
      }
      if (!filterForTOFEventTime(t)) {
        continue;
      }
      histos.fill(HIST("goodforevtime/p"), t.p());
      histos.fill(HIST("goodforevtime/pt"), t.pt());
      histos.fill(HIST("goodforevtime/length"), t.length());
      histos.fill(HIST("goodforevtime/tofSignal"), t.tofSignal());
      histos.fill(HIST("goodforevtime/beta"), t.p(), beta);
      histos.fill(HIST("goodforevtime/delta"), t.p(), t.tofSignal() - t.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(t));
      histos.fill(HIST("goodforevtime/expP"), t.p(), t.tofExpMom());
      histos.fill(HIST("goodforevtime/mass"), mass);
      histos.fill(HIST("goodforevtime/tofSignalPerCollision"), ncolls % 6000, t.tofSignal());
    }
  }
};

/// Task to produce the TOF QA plots
struct tofPidFullQa {
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

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

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

  template <uint8_t i>
  void addParticleHistos(const AxisSpec& pAxis, const AxisSpec& ptAxis)
  {
    // Exp signal
    const AxisSpec expAxis{1000, 0, 2e6, Form("t_{exp}(%s) (ps)", pT[i])};
    histos.add(hexpected[i].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("t-t_{ev}-t_{exp}(%s) (ps)", pT[i])};
    histos.add(hexpected_diff[i].data(), "", kTH2F, {pAxis, deltaAxis});

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TOF}(%s) (ps)", pT[i])};
    histos.add(hexpsigma[i].data(), "", kTH2F, {pAxis, expSigmaAxis});

    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TOF}(%s)", pT[i]);
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
    histos.add(hnsigma[i].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigmapt[i].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmapospt[i].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmanegpt[i].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
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

    h = histos.add<TH1>("event/trackselection", "", kTH1F, {{10, 0, 10, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(3, "hasITS");
    h->GetXaxis()->SetBinLabel(4, "hasTOF");

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
      addParticleHistos<i>(pAxis, ptAxis);
    });
  }

  template <o2::track::PID::ID id, typename T>
  void fillParticleHistos(const T& t, const float& tof)
  {
    if (applyRapidityCut) {
      const float y = TMath::ASinH(t.pt() / TMath::Sqrt(PID::getMass2(id) + t.pt() * t.pt()) * TMath::SinH(t.eta()));
      if (abs(y) > 0.5) {
        return;
      }
    }

    const auto& nsigma = o2::aod::pidutils::tofNSigma<id>(t);
    const auto& diff = o2::aod::pidutils::tofExpSignalDiff<id>(t);
    histos.fill(HIST(hexpected[id]), t.p(), tof - diff);
    histos.fill(HIST(hexpected_diff[id]), t.p(), diff);
    histos.fill(HIST(hexpsigma[id]), t.p(), o2::aod::pidutils::tofExpSigma<id>(t));
    histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
    histos.fill(HIST(hnsigmapt[id]), t.pt(), nsigma);
    if (t.sign() > 0) {
      histos.fill(HIST(hnsigmapospt[id]), t.pt(), nsigma);
    } else {
      histos.fill(HIST(hnsigmanegpt[id]), t.pt(), nsigma);
    }
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                         aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                         aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                         aod::TOFSignal, aod::TrackSelection>;
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               Trks const& tracks)
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
    int tofmult = 0;
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
    // if (0 && ntracks < 1) {
    //   return;
    // }
    // if (0 && tofmult < 1) {
    //   return;
    // }
    histos.fill(HIST("event/evsel"), 3);
    if (abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("event/evsel"), 4);
    histos.fill(HIST("event/vertexz"), collision.posZ());
    histos.fill(HIST("event/trackmultiplicity"), ntracks);
    histos.fill(HIST("event/tofmultiplicity"), tofmult);

    const float collisionTime_ps = collision.collisionTime() * 1000.f;
    histos.fill(HIST("event/colltime"), collisionTime_ps);
    histos.fill(HIST("event/colltimereso"), tofmult, collision.collisionTimeRes() * 1000.f);

    for (auto t : tracks) {
      histos.fill(HIST("event/trackselection"), 0.5f);
      if (!t.isGlobalTrack()) { // Skipping non global tracks
        continue;
      }
      histos.fill(HIST("event/trackselection"), 1.5f);
      if (!t.hasITS()) { // Skipping tracks without ITS
        continue;
      }
      histos.fill(HIST("event/trackselection"), 2.5f);
      if (!t.hasTOF()) { // Skipping tracks without TOF
        continue;
      }
      histos.fill(HIST("event/trackselection"), 3.5f);

      const float tof = t.tofSignal() - collisionTime_ps;

      histos.fill(HIST("event/particlehypo"), t.pidForTracking());
      histos.fill(HIST("event/tofsignal"), t.p(), t.tofSignal());
      histos.fill(HIST("event/pexp"), t.p(), t.tofExpMom());
      histos.fill(HIST("event/eta"), t.eta());
      histos.fill(HIST("event/phi"), t.phi());
      histos.fill(HIST("event/etaphi"), t.eta(), t.phi());
      histos.fill(HIST("event/length"), t.length());
      histos.fill(HIST("event/pt"), t.pt());
      histos.fill(HIST("event/p"), t.p());
      // histos.fill(HIST("event/ptreso"), t.p(), t.sigma1Pt() * t.pt() * t.pt());
      //
      fillParticleHistos<PID::Electron>(t, tof);
      fillParticleHistos<PID::Muon>(t, tof);
      fillParticleHistos<PID::Pion>(t, tof);
      fillParticleHistos<PID::Kaon>(t, tof);
      fillParticleHistos<PID::Proton>(t, tof);
      fillParticleHistos<PID::Deuteron>(t, tof);
      fillParticleHistos<PID::Triton>(t, tof);
      fillParticleHistos<PID::Helium3>(t, tof);
      fillParticleHistos<PID::Alpha>(t, tof);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofSignal>(cfgc),
                               adaptAnalysisTask<tofPidFull>(cfgc)};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tofPidFullQa>(cfgc));
  }
  if (cfgc.options().get<int>("add-qa-ev-time")) {
    workflow.push_back(adaptAnalysisTask<tofPidCollisionTimeQa>(cfgc));
  }
  return workflow;
}
