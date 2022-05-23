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
/// \file   pidTOFBase.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TOF PID tasks.
///

#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include <CCDB/BasicCCDBManager.h>
#include "TableHelper.h"
#include "TOFBase/EventTimeMaker.h"
#include "pidTOFBase.h"

using namespace o2;
using namespace o2::pid;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TOF PID QA histograms for TOF event time"}},
                                       {"evtime", VariantType::Int, 1, {"Produce the table for the Event Time"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the TOF signal from the trackTime information
struct tofSignal {
  o2::framework::Produces<o2::aod::TOFSignal> table;
  bool enableTable = false;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "TOFSignal");
    if (enableTable) {
      LOG(info) << "Table TOFSignal enabled!";
    }
  }
  using Trks = o2::soa::Join<aod::Tracks, aod::TracksExtra>;
  void process(Trks const& tracks)
  {
    if (!enableTable) {
      return;
    }
    table.reserve(tracks.size());
    for (auto& t : tracks) {
      table(o2::pid::tof::TOFSignal<Trks::iterator>::GetTOFSignal(t));
    }
  }
};

/// Selection criteria for tracks used for TOF event time
template <typename trackType>
bool filterForTOFEventTime(const trackType& tr)
{
  return (tr.hasTOF() && tr.p() > 0.5f && tr.p() < 2.f && tr.trackType() == o2::aod::track::TrackTypeEnum::Track);
} // accept all

/// Specialization of TOF event time maker
template <typename trackType,
          bool (*trackFilter)(const trackType&),
          template <typename T, o2::track::PID::ID> typename response,
          typename trackTypeContainer,
          typename responseParametersType>
o2::tof::eventTimeContainer evTimeMakerForTracks(const trackTypeContainer& tracks,
                                                 const responseParametersType& responseParameters,
                                                 const float& diamond = 6.0)
{
  return o2::tof::evTimeMakerFromParam<trackTypeContainer, trackType, trackFilter, response, responseParametersType>(tracks, responseParameters, diamond);
}

/// Task to produce the TOF event time table
struct tofEventTime {
  // Tables to produce
  Produces<o2::aod::TOFEvTime> tableEvTime;
  Produces<o2::aod::pidEvTimeFlags> tableFlags;
  static constexpr bool removeTOFEvTimeBias = true; // Flag to substract the Ev. Time bias for low multiplicity events with TOF
  static constexpr float diamond = 6.0;             // Collision diamond used in the estimation of the TOF event time

  bool enableTable = false;
  // Detector response and input parameters
  DetectorResponse response;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> sigmaname{"param-sigma", "TOFReso", "Name of the parametrization for the expected sigma, used in both file and CCDB mode"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TOF", "Path of the TOF parametrization on the CCDB"};
  Configurable<long> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};

  void init(o2::framework::InitContext& initContext)
  {
    // Check that both processes are not enabled
    if (doprocessNoFT0 == true && doprocessFT0 == true) {
      LOGF(fatal, "Cannot enable processNoEvTime and processEvTime at the same time. Please choose one.");
    }
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "TOFEvTime");
    if (!enableTable) {
      LOG(info) << "Table for TOF Event time (TOFEvTime) is not required, disabling it";
      return;
    }
    LOG(info) << "Table TOFEvTime enabled!";

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

  using TrksEvTime = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal>;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksEvTime::iterator, pid>;
  void processNoFT0(TrksEvTime const& tracks,
                    aod::Collisions const&)
  {
    if (!enableTable) {
      return;
    }

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());

    int lastCollisionId = -1;      // Last collision ID analysed
    for (auto const& t : tracks) { // Loop on collisions
      if (!t.has_collision()) {    // Track was not assigned, cannot compute event time
        tableFlags(0);
        tableEvTime(0.f, 999.f, -1);
        continue;
      }
      if (t.collisionId() == lastCollisionId) { // Event time from this collision is already in the table
        continue;
      }
      /// Create new table for the tracks in a collision
      lastCollisionId = t.collisionId(); /// Cache last collision ID

      const auto& tracksInCollision = tracks.sliceBy(aod::track::collisionId, lastCollisionId);
      // First make table for event time
      const auto evTimeTOF = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, response, diamond);
      int nGoodTracksForTOF = 0;
      float et = evTimeTOF.mEventTime;
      float erret = evTimeTOF.mEventTimeError;

      for (auto const& trk : tracksInCollision) { // Loop on Tracks
        if constexpr (removeTOFEvTimeBias) {
          evTimeTOF.removeBias<TrksEvTime::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, et, erret, 2);
        }
        uint8_t flags = 0;
        if (erret > 199.f) {
          flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeTOF;
        }
        tableFlags(flags);
        tableEvTime(et, erret, evTimeTOF.mEventTimeMultiplicity);
      }
    }
  }

  PROCESS_SWITCH(tofEventTime, processNoFT0, "Process without FT0", true);

  using EvTimeCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>;
  void processFT0(TrksEvTime const& tracks,
                  aod::FT0s const&,
                  EvTimeCollisions const&)
  {
    if (!enableTable) {
      return;
    }

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());

    int lastCollisionId = -1;      // Last collision ID analysed
    for (auto const& t : tracks) { // Loop on collisions
      if (!t.has_collision()) {    // Track was not assigned, cannot compute event time
        tableFlags(0);
        tableEvTime(0.f, 999.f, -1);
        continue;
      }
      if (t.collisionId() == lastCollisionId) { // Event time from this collision is already in the table
        continue;
      }
      /// Create new table for the tracks in a collision
      lastCollisionId = t.collisionId(); /// Cache last collision ID

      const auto& tracksInCollision = tracks.sliceBy(aod::track::collisionId, lastCollisionId);
      const auto& collision = t.collision_as<EvTimeCollisions>();

      // Compute the TOF event time
      const auto evTimeTOF = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, response, diamond);

      float t0AC[2] = {.0f, 999.f};                                       // Value and error of T0A or T0C or T0AC
      float t0TOF[2] = {evTimeTOF.mEventTime, evTimeTOF.mEventTimeError}; // Value and error of TOF

      uint8_t flags = 0;
      int nGoodTracksForTOF = 0;
      float eventTime = 0.f;
      float sumOfWeights = 0.f;
      float weight = 0.f;
      for (auto const& trk : tracksInCollision) { // Loop on Tracks
        // Reset the flag
        flags = 0;
        // Reset the event time
        eventTime = 0.f;
        sumOfWeights = 0.f;
        weight = 0.f;
        // Remove the bias on TOF ev. time
        if constexpr (removeTOFEvTimeBias) {
          evTimeTOF.removeBias<TrksEvTime::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, t0TOF[0], t0TOF[1], 2);
        }
        if (t0TOF[1] < 199.f) {
          flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeTOF;

          weight = 1.f / (t0TOF[1] * t0TOF[1]);
          eventTime += t0TOF[0] * weight;
          sumOfWeights += weight;
        }

        if (collision.has_foundFT0()) { // T0 measurement is available
          // const auto& ft0 = collision.foundFT0();
          if (collision.t0ACValid()) {
            t0AC[0] = collision.t0AC();
            t0AC[1] = collision.t0resolution();
            flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeT0AC;
          }

          weight = 1.f / (t0AC[1] * t0AC[1]);
          eventTime += t0AC[0] * weight;
          sumOfWeights += weight;
        }

        tableFlags(flags);
        tableEvTime(eventTime / sumOfWeights, sqrt(1. / sumOfWeights), evTimeTOF.mEventTimeMultiplicity);
      }
    }
  }

  PROCESS_SWITCH(tofEventTime, processFT0, "Process with FT0", false);
};

/// Task that checks the TOF collision time
struct tofPidCollisionTimeQa {
  Configurable<int> nBinsEvTime{"nBinsEvTime", 1000, "Number of bins for the event time"};
  Configurable<float> minEvTime{"minEvTime", -1000.f, "Minimum in range in event time"};
  Configurable<float> maxEvTime{"maxEvTime", 1000.f, "Maximum in range in event time"};
  Configurable<int> nBinsTofSignal{"nBinsTofSignal", 5000, "Number of bins for the tof signal time"};
  Configurable<float> minTofSignal{"minTofSignal", 0.f, "Minimum in range in tof signal time"};
  Configurable<float> maxTofSignal{"maxTofSignal", 100e3, "Maximum in range in tof signal time"};
  Configurable<int> nBinsEvTimeReso{"nBinsEvTimeReso", 1000, "Number of bins in event time resolution"};
  Configurable<float> rangeEvTimeReso{"rangeEvTimeReso", 1000.f, "Range in event time resolution"};
  Configurable<int> nBinsMultiplicity{"nBinsMultiplicity", 1000, "Number of bins for the multiplicity"};
  Configurable<float> rangeMultiplicity{"rangeMultiplicity", 1000.f, "Range for the multiplicity"};
  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 200, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext& initContext)
  {
    const AxisSpec evTimeAxis{nBinsEvTime, minEvTime, maxEvTime, "TOF event time (ps)"};
    const AxisSpec multAxis{nBinsEvTime, 0, rangeMultiplicity, "Track multiplicity for TOF event time"};
    const AxisSpec evTimeResoAxis{nBinsEvTimeReso, 0, rangeEvTimeReso, "TOF event time resolution (ps)"};
    const AxisSpec tofSignalAxis{nBinsTofSignal, minTofSignal, maxTofSignal, "TOF signal (ps)"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} GeV/#it{c}"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} GeV/#it{c}"};
    if (logAxis) {
      pAxis.makeLogaritmic();
      ptAxis.makeLogaritmic();
    }
    const AxisSpec collisionAxis{6000, -0.5f, 6000.f - .5f, "Collision index % 6000"};
    const AxisSpec massAxis{1000, 0, 3, "TOF mass (GeV/#it{c}^{2})"};
    const AxisSpec betaAxis{1000, 0, 1.5, "TOF #beta"};
    const AxisSpec deltaAxis{1000, -10000, 10000, "t-t_{ev}-t_{exp}(#pi) (ps)"};
    const AxisSpec lengthAxis{1000, 0, 600, "Track length (cm)"};

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
    histos.add("eventTime", "eventTime", kTH1F, {evTimeAxis});
    histos.add("eventTimeReso", "eventTimeReso", kTH1F, {evTimeResoAxis});
    histos.add("eventTimeMult", "eventTimeMult", kTH1F, {multAxis});
    histos.add("eventTimeVsMult", "eventTimeVsMult", kTH2F, {multAxis, evTimeAxis});
    histos.add("eventTimeResoVsMult", "eventTimeResoVsMult", kTH2F, {multAxis, evTimeResoAxis});
    histos.add<TH1>("collisionTime", "collisionTime", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("Collision time (ps)");
    histos.add<TH1>("collisionTimeRes", "collisionTimeRes", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("Collision time resolution (ps)");

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

    histos.add("goodreso/p", "p", kTH1F, {pAxis});
    histos.add("goodreso/pt", "pt", kTH1F, {ptAxis});
    histos.add("goodreso/ptden", "ptden", kTH1F, {ptAxis});
    histos.add("goodreso/length", "length", kTH1F, {lengthAxis});
    histos.add("goodreso/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("goodreso/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("goodreso/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("goodreso/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("goodreso/mass", "mass", kTH1F, {massAxis});
    histos.add("goodreso/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("badreso/p", "p", kTH1F, {pAxis});
    histos.add("badreso/pt", "pt", kTH1F, {ptAxis});
    histos.add("badreso/ptden", "ptden", kTH1F, {ptAxis});
    histos.add("badreso/length", "length", kTH1F, {lengthAxis});
    histos.add("badreso/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("badreso/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("badreso/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("badreso/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("badreso/mass", "mass", kTH1F, {massAxis});
    histos.add("badreso/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("goodforevtime/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("goodforevtime/p", "p", kTH1F, {pAxis});
    histos.add("goodforevtime/pt", "pt", kTH1F, {ptAxis});
    histos.add("goodforevtime/length", "length", kTH1F, {lengthAxis});
    histos.add("goodforevtime/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("goodforevtime/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("goodforevtime/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("goodforevtime/mass", "mass", kTH1F, {massAxis});
    histos.add("goodforevtime/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("withqualitycuts/p", "p", kTH1F, {pAxis});
    histos.add("withqualitycuts/pt", "pt", kTH1F, {ptAxis});
    histos.add("withqualitycuts/length", "length", kTH1F, {lengthAxis});
    histos.add("withqualitycuts/mass", "mass", kTH1F, {massAxis});
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection>;
  void process(Trks const& tracks, aod::Collisions const&)
  {
    static int ncolls = 0;
    int lastCollisionId = -1; // Last collision ID analysed
    for (auto& t : tracks) {
      if (!t.has_collision()) { // Track was not assigned to a collision
        continue;
      } else if (t.collisionId() == lastCollisionId) { // Event was already processed
        continue;
      }
      /// Create new table for the tracks in a collision
      lastCollisionId = t.collisionId(); /// Cache last collision ID

      // const auto& collision = t.collision();
      histos.fill(HIST("eventSelection"), 0.5f);
      histos.fill(HIST("eventSelection"), 1.5f);
      if (t.tofEvTimeErr() > 199.f) {
        histos.fill(HIST("eventSelection"), 2.5f);
      } else {
        histos.fill(HIST("eventSelection"), 3.5f);
      }
      histos.fill(HIST("eventTime"), t.tofEvTime());
      histos.fill(HIST("eventTimeReso"), t.tofEvTimeErr());
      histos.fill(HIST("eventTimeMult"), t.tofEvTimeMult());
      histos.fill(HIST("eventTimeVsMult"), t.tofEvTimeMult(), t.tofEvTime());
      histos.fill(HIST("eventTimeResoVsMult"), t.tofEvTimeMult(), t.tofEvTimeErr());

      histos.fill(HIST("collisionTime"), t.collision().collisionTime());
      histos.fill(HIST("collisionTimeRes"), t.collision().collisionTimeRes());
      ncolls++;

      const auto tracksInCollision = tracks.sliceBy(aod::track::collisionId, lastCollisionId);

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

        if (trk.tofEvTimeErr() > 199.f) {
          histos.fill(HIST("badreso/ptden"), trk.pt());
        } else {
          histos.fill(HIST("goodreso/ptden"), trk.pt());
        }

        if (!trk.hasTOF()) {
          continue;
        }
        histos.fill(HIST("trackSelection"), 4.5f);

        const float beta = o2::pid::tof::Beta<Trks::iterator>::GetBeta(trk, trk.tofEvTime());
        const float mass = o2::pid::tof::TOFMass<Trks::iterator>::GetTOFMass(trk.p(), beta);
        histos.fill(HIST("withtof/p"), trk.p());
        histos.fill(HIST("withtof/pt"), trk.pt());
        histos.fill(HIST("withtof/length"), trk.length());
        histos.fill(HIST("withtof/tofSignal"), trk.tofSignal());
        histos.fill(HIST("withtof/beta"), trk.p(), beta);
        histos.fill(HIST("withtof/delta"), trk.p(), trk.tofSignal() - trk.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(trk));
        histos.fill(HIST("withtof/expP"), trk.p(), trk.tofExpMom());
        histos.fill(HIST("withtof/mass"), mass);
        histos.fill(HIST("withtof/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
        if (trk.pt() > 0.3 && beta > 0.3) {
          histos.fill(HIST("withqualitycuts/p"), trk.p());
          histos.fill(HIST("withqualitycuts/pt"), trk.pt());
          histos.fill(HIST("withqualitycuts/length"), trk.length());
          histos.fill(HIST("withqualitycuts/mass"), mass);
        }

        if (trk.tofEvTimeErr() > 199.f) {
          histos.fill(HIST("badreso/p"), trk.p());
          histos.fill(HIST("badreso/pt"), trk.pt());
          histos.fill(HIST("badreso/length"), trk.length());
          histos.fill(HIST("badreso/tofSignal"), trk.tofSignal());
          histos.fill(HIST("badreso/beta"), trk.p(), beta);
          histos.fill(HIST("badreso/delta"), trk.p(), trk.tofSignal() - trk.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(trk));
          histos.fill(HIST("badreso/expP"), trk.p(), trk.tofExpMom());
          histos.fill(HIST("badreso/mass"), mass);
          histos.fill(HIST("badreso/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
        } else {
          histos.fill(HIST("goodreso/p"), trk.p());
          histos.fill(HIST("goodreso/pt"), trk.pt());
          histos.fill(HIST("goodreso/length"), trk.length());
          histos.fill(HIST("goodreso/tofSignal"), trk.tofSignal());
          histos.fill(HIST("goodreso/beta"), trk.p(), beta);
          histos.fill(HIST("goodreso/delta"), trk.p(), trk.tofSignal() - trk.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(trk));
          histos.fill(HIST("goodreso/expP"), trk.p(), trk.tofExpMom());
          histos.fill(HIST("goodreso/mass"), mass);
          histos.fill(HIST("goodreso/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
        }
        if (!filterForTOFEventTime(trk)) {
          continue;
        }
        histos.fill(HIST("goodforevtime/p"), trk.p());
        histos.fill(HIST("goodforevtime/pt"), trk.pt());
        histos.fill(HIST("goodforevtime/length"), trk.length());
        histos.fill(HIST("goodforevtime/tofSignal"), trk.tofSignal());
        histos.fill(HIST("goodforevtime/beta"), trk.p(), beta);
        histos.fill(HIST("goodforevtime/delta"), trk.p(), trk.tofSignal() - trk.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(trk));
        histos.fill(HIST("goodforevtime/expP"), trk.p(), trk.tofExpMom());
        histos.fill(HIST("goodforevtime/mass"), mass);
        histos.fill(HIST("goodforevtime/tofSignalPerCollision"), ncolls % 6000, trk.tofSignal());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofSignal>(cfgc)};
  if (!cfgc.options().get<int>("evtime")) {
    return workflow;
  }
  workflow.push_back(adaptAnalysisTask<tofEventTime>(cfgc));
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tofPidCollisionTimeQa>(cfgc));
  }
  return workflow;
}
