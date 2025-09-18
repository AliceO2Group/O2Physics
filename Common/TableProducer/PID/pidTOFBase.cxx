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

#include "pidTOFBase.h"

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseTOF.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsTOF/ParameterContainers.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <PID/PIDTOF.h>
#include <ReconstructionDataFormats/PID.h>
#include <TOFBase/EventTimeMaker.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Selection criteria for tracks used for TOF event time
float trackDistanceForGoodMatch = 999.f;
float trackDistanceForGoodMatchLowMult = 999.f;
int multiplicityThreshold = 0;
using Run3Trks = o2::soa::Join<aod::TracksIU, aod::TracksExtra>;
using Run3Cols = aod::Collisions;
bool isTrackGoodMatchForTOFPID(const Run3Trks::iterator& tr, const Run3Cols& /*ev*/)
{
  if (!tr.hasTOF()) {
    return false;
  }
  return true;
}

/// Task to produce the TOF signal from the trackTime information
struct tofSignal {
  o2::framework::Produces<o2::aod::TOFSignal> table;
  o2::framework::Produces<o2::aod::pidTOFFlags> tableFlags;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  bool enableTable = false;      // Flag to check if the TOF signal table is requested or not
  bool enableTableFlags = false; // Flag to check if the TOF signal flags table is requested or not
  // CCDB configuration
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"};
  Configurable<float> distanceForGoodMatch{"distanceForGoodMatch", 999.f, "Maximum distance to consider a good match"};
  Configurable<float> distanceForGoodMatchLowMult{"distanceForGoodMatchLowMult", 999.f, "Maximum distance to consider a good match for low multiplicity events"};
  Configurable<int> multThreshold{"multThreshold", 0, "Multiplicity threshold to consider a low multiplicity event"};
  Configurable<bool> enableQaHistograms{"enableQaHistograms", false, "Flag to enable the QA histograms"};

  void init(o2::framework::InitContext& initContext)
  {
    if (doprocessRun2 && doprocessRun3) {
      LOG(fatal) << "Both processRun2 and processRun3 are enabled. Pick one of the two";
    }
    if (!doprocessRun2 && !doprocessRun3) {
      LOG(fatal) << "Neither processRun2 nor processRun3 are enabled. Pick one of the two";
    }

    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "TOFSignal");
    if (enableTable) {
      LOG(info) << "Table TOFSignal enabled!";
    }
    enableTableFlags = isTableRequiredInWorkflow(initContext, "pidTOFFlags");
    if (enableTableFlags) {
      LOG(info) << "Table pidTOFFlags enabled!";
    }
    trackDistanceForGoodMatch = distanceForGoodMatch;
    trackDistanceForGoodMatchLowMult = distanceForGoodMatchLowMult;
    multiplicityThreshold = multThreshold;
    LOG(info) << "Configuring selections for good match: " << trackDistanceForGoodMatch << " low mult " << trackDistanceForGoodMatchLowMult << " mult. threshold " << multiplicityThreshold;
    if (!enableQaHistograms) {
      return;
    }
    histos.add("tofSignal", "tofSignal", kTH1D, {{1000, -1000, 1000000, "tofSignal (ps)"}});
    if (enableTableFlags) {
      histos.add("goodForPIDFlags", "goodForPIDFlags", kTH1D, {{3, 0, 3, "flags"}});
    }
  }
  void processRun3(Run3Trks const& tracks, Run3Cols const& collisions)
  {
    if (!enableTable) {
      return;
    }
    table.reserve(tracks.size());
    if (enableTableFlags) {
      tableFlags.reserve(tracks.size());
    }
    for (const auto& t : tracks) {
      const auto s = o2::pid::tof::TOFSignal<Run3Trks::iterator>::GetTOFSignal(t);
      if (enableQaHistograms) {
        histos.fill(HIST("tofSignal"), s);
      }
      table(s);
      if (!enableTableFlags) {
        continue;
      }
      const auto b = isTrackGoodMatchForTOFPID(t, collisions);
      if (enableQaHistograms) {
        histos.fill(HIST("goodForPIDFlags"), s);
      }
      tableFlags(b);
    }
  }
  PROCESS_SWITCH(tofSignal, processRun3, "Process Run3 data i.e. input is TrackIU", true);

  using TrksRun2 = o2::soa::Join<aod::Tracks, aod::TracksExtra>;
  void processRun2(TrksRun2 const& tracks)
  {
    if (!enableTable) {
      return;
    }
    table.reserve(tracks.size());
    if (enableTableFlags) {
      tableFlags.reserve(tracks.size());
    }
    for (const auto& t : tracks) {
      table(o2::pid::tof::TOFSignal<TrksRun2::iterator>::GetTOFSignal(t));
      if (!enableTableFlags) {
        continue;
      }
      tableFlags(true);
    }
  }
  PROCESS_SWITCH(tofSignal, processRun2, "Process Run2 data i.e. input is Tracks", false);
};

/// Selection criteria for tracks used for TOF event time
float trackSampleMinMomentum = 0.5f;
float trackSampleMaxMomentum = 2.f;
template <typename trackType>
bool filterForTOFEventTime(const trackType& tr)
{
  return (tr.hasTOF() && tr.p() > trackSampleMinMomentum && tr.p() < trackSampleMaxMomentum && (tr.trackType() == o2::aod::track::TrackTypeEnum::Track || tr.trackType() == o2::aod::track::TrackTypeEnum::TrackIU));
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
  Produces<o2::aod::EvTimeTOFOnly> tableEvTimeTOFOnly;
  Produces<o2::aod::pidEvTimeFlags> tableFlags;
  static constexpr bool removeTOFEvTimeBias = true; // Flag to subtract the Ev. Time bias for low multiplicity events with TOF
  static constexpr float diamond = 6.0;             // Collision diamond used in the estimation of the TOF event time
  static constexpr float errDiamond = diamond * 33.356409f;
  static constexpr float weightDiamond = 1.f / (errDiamond * errDiamond);

  bool enableTable = false;
  bool enableTableTOFOnly = false;
  // Detector response and input parameters
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<bool> inheritFromBaseTask{"inheritFromBaseTask", true, "Flag to iherit all common configurables from the TOF base task"};
  // CCDB configuration (inherited from TOF signal task)
  Configurable<std::string> url{"ccdb-url", "", "url of the ccdb repository"};
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  // Event time configurations
  Configurable<float> minMomentum{"minMomentum", 0.5f, "Minimum momentum to select track sample for TOF event time"};
  Configurable<float> maxMomentum{"maxMomentum", 2.0f, "Maximum momentum to select track sample for TOF event time"};
  Configurable<float> maxEvTimeTOF{"maxEvTimeTOF", 100000.0f, "Maximum value of the TOF event time"};
  Configurable<bool> sel8TOFEvTime{"sel8TOFEvTime", false, "Flag to compute the ev. time only for events that pass the sel8 ev. selection"};
  Configurable<int> maxNtracksInSet{"maxNtracksInSet", 10, "Size of the set to consider for the TOF ev. time computation"};
  // TOF Calib configuration
  Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
  Configurable<std::string> parametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<bool> loadResponseFromCCDB{"loadResponseFromCCDB", false, "Flag to load the response from the CCDB"};
  Configurable<bool> enableTimeDependentResponse{"enableTimeDependentResponse", false, "Flag to use the collision timestamp to fetch the PID Response"};
  Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};

  void init(o2::framework::InitContext& initContext)
  {
    if (inheritFromBaseTask.value) {
      if (!getTaskOptionValue(initContext, "tof-signal", "ccdb-url", url.value, true)) {
        LOG(fatal) << "Could not get ccdb-url from tof-signal task";
      }
      if (!getTaskOptionValue(initContext, "tof-signal", "ccdb-timestamp", timestamp.value, true)) {
        LOG(fatal) << "Could not get ccdb-timestamp from tof-signal task";
      }
    }

    trackSampleMinMomentum = minMomentum;
    trackSampleMaxMomentum = maxMomentum;
    LOG(info) << "Configuring track sample for TOF ev. time: " << trackSampleMinMomentum << " < p < " << trackSampleMaxMomentum;
    // Check that both processes are not enabled
    int nEnabled = 0;
    if (doprocessRun2 == true) {
      LOGF(info, "Enabling process function: processRun2");
      nEnabled++;
    }
    if (doprocessNoFT0 == true) {
      LOGF(info, "Enabling process function: processNoFT0");
      nEnabled++;
    }
    if (doprocessFT0 == true) {
      LOGF(info, "Enabling process function: processFT0");
      nEnabled++;
    }
    if (doprocessOnlyFT0 == true) {
      LOGF(info, "Enabling process function: processOnlyFT0");
      nEnabled++;
    }
    if (nEnabled > 1) {
      LOGF(fatal, "Cannot enable more process functions at the same time. Please choose one.");
    }
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "TOFEvTime");
    if (!enableTable) {
      LOG(info) << "Table for TOF Event time (TOFEvTime) is not required, disabling it";
      return;
    }
    LOG(info) << "Table TOFEvTime enabled!";

    enableTableTOFOnly = isTableRequiredInWorkflow(initContext, "EvTimeTOFOnly");
    if (enableTableTOFOnly) {
      LOG(info) << "Table EvTimeTOFOnly enabled!";
    }

    if (sel8TOFEvTime.value == true) {
      LOG(info) << "TOF event time will be computed for collisions that pass the event selection only!";
    }
    // Getting the parametrization parameters
    ccdb->setURL(url.value);
    ccdb->setTimestamp(timestamp.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    //

    // TODO: implement the automatic pass name detection from metadata
    if (passName.value == "") {
      passName.value = "unanchored"; // temporary default
      LOG(warning) << "Passed autodetect mode for pass, not implemented yet, waiting for metadata. Taking '" << passName.value << "'";
    }
    LOG(info) << "Using parameter collection, starting from pass '" << passName.value << "'";

    const std::string fname = paramFileName.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file " << fname << ", using param: " << parametrizationPath.value;
      if (1) {
        o2::tof::ParameterCollection paramCollection;
        paramCollection.loadParamFromFile(fname, parametrizationPath.value);
        LOG(info) << "+++ Loaded parameter collection from file +++";
        if (!paramCollection.retrieveParameters(mRespParamsV2, passName.value)) {
          if (fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          }
        } else {
          mRespParamsV2.setShiftParameters(paramCollection.getPars(passName.value));
          mRespParamsV2.printShiftParameters();
        }
      } else {
        mRespParamsV2.loadParamFromFile(fname.data(), parametrizationPath.value);
      }
    } else if (loadResponseFromCCDB) { // Loading it from CCDB
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << parametrizationPath.value << " for timestamp " << timestamp.value;

      o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, timestamp.value);
      paramCollection->print();
      if (!paramCollection->retrieveParameters(mRespParamsV2, passName.value)) {
        if (fatalOnPassNotAvailable) {
          LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        } else {
          LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        }
      } else {
        mRespParamsV2.setShiftParameters(paramCollection->getPars(passName.value));
        mRespParamsV2.printShiftParameters();
      }
    }
    mRespParamsV2.print();
    o2::tof::eventTimeContainer::setMaxNtracksInSet(maxNtracksInSet.value);
    o2::tof::eventTimeContainer::printConfig();
  }

  ///
  /// Process function to prepare the event for each track on Run 2 data
  void processRun2(aod::Tracks const& tracks,
                   aod::Collisions const&)
  {
    if (!enableTable) {
      return;
    }

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());

    for (auto const& t : tracks) { // Loop on collisions
      if (!t.has_collision()) {    // Track was not assigned, cannot compute event time
        tableFlags(0);
        tableEvTime(0.f, 999.f);
        continue;
      }
      tableFlags(1);
      tableEvTime(t.collision().collisionTime() * 1000.f, t.collision().collisionTimeRes() * 1000.f);
    }
  }
  PROCESS_SWITCH(tofEventTime, processRun2, "Process with Run2 data", false);

  ///
  /// Process function to prepare the event for each track on Run 3 data without the FT0
  using TrksEvTime = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TOFSignal>;
  // Define slice per collision
  Preslice<TrksEvTime> perCollision = aod::track::collisionId;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksEvTime::iterator, pid>;
  using EvTimeCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  void processNoFT0(TrksEvTime const& tracks,
                    EvTimeCollisions const&)
  {
    if (!enableTable) {
      return;
    }

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());
    if (enableTableTOFOnly) {
      tableEvTimeTOFOnly.reserve(tracks.size());
    }

    int lastCollisionId = -1;                                                                                    // Last collision ID analysed
    for (auto const& t : tracks) {                                                                               // Loop on collisions
      if (!t.has_collision() || ((sel8TOFEvTime.value == true) && !t.collision_as<EvTimeCollisions>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
        tableFlags(0);
        tableEvTime(0.f, 999.f);
        if (enableTableTOFOnly) {
          tableEvTimeTOFOnly((uint8_t)0, 0.f, 0.f, -1);
        }
        continue;
      }
      if (t.collisionId() == lastCollisionId) { // Event time from this collision is already in the table
        continue;
      }
      /// Create new table for the tracks in a collision
      lastCollisionId = t.collisionId(); /// Cache last collision ID

      const auto& tracksInCollision = tracks.sliceBy(perCollision, lastCollisionId);

      // First make table for event time
      const auto evTimeTOF = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV2, diamond);
      int nGoodTracksForTOF = 0;
      float et = evTimeTOF.mEventTime;
      float erret = evTimeTOF.mEventTimeError;

      for (auto const& trk : tracksInCollision) { // Loop on Tracks
        if constexpr (removeTOFEvTimeBias) {
          evTimeTOF.removeBias<TrksEvTime::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, et, erret, 2);
        }
        uint8_t flags = 0;
        if (erret < errDiamond && (maxEvTimeTOF <= 0.f || std::abs(et) < maxEvTimeTOF)) {
          flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeTOF;
        } else {
          et = 0.f;
          erret = errDiamond;
        }
        tableFlags(flags);
        tableEvTime(et, erret);
        if (enableTableTOFOnly) {
          tableEvTimeTOFOnly((uint8_t)filterForTOFEventTime(trk), et, erret, evTimeTOF.mEventTimeMultiplicity);
        }
      }
    }
  }
  PROCESS_SWITCH(tofEventTime, processNoFT0, "Process without FT0", true);

  ///
  /// Process function to prepare the event for each track on Run 3 data with the FT0
  using EvTimeCollisionsFT0 = soa::Join<EvTimeCollisions, aod::FT0sCorrected>;
  void processFT0(TrksEvTime const& tracks,
                  aod::FT0s const&,
                  EvTimeCollisionsFT0 const&)
  {
    if (!enableTable) {
      return;
    }

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());
    if (enableTableTOFOnly) {
      tableEvTimeTOFOnly.reserve(tracks.size());
    }

    int lastCollisionId = -1;                                                                                       // Last collision ID analysed
    for (auto const& t : tracks) {                                                                                  // Loop on collisions
      if (!t.has_collision() || ((sel8TOFEvTime.value == true) && !t.collision_as<EvTimeCollisionsFT0>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
        tableFlags(0);
        tableEvTime(0.f, 999.f);
        if (enableTableTOFOnly) {
          tableEvTimeTOFOnly((uint8_t)0, 0.f, 0.f, -1);
        }
        continue;
      }
      if (t.collisionId() == lastCollisionId) { // Event time from this collision is already in the table
        continue;
      }
      /// Create new table for the tracks in a collision
      lastCollisionId = t.collisionId(); /// Cache last collision ID

      const auto& tracksInCollision = tracks.sliceBy(perCollision, lastCollisionId);
      const auto& collision = t.collision_as<EvTimeCollisionsFT0>();

      // Compute the TOF event time
      const auto evTimeTOF = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV2, diamond);

      float t0AC[2] = {.0f, 999.f};                                                                                   // Value and error of T0A or T0C or T0AC
      float t0TOF[2] = {static_cast<float_t>(evTimeTOF.mEventTime), static_cast<float_t>(evTimeTOF.mEventTimeError)}; // Value and error of TOF

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
        if (t0TOF[1] < errDiamond && (maxEvTimeTOF <= 0 || std::abs(t0TOF[0]) < maxEvTimeTOF)) {
          flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeTOF;

          weight = 1.f / (t0TOF[1] * t0TOF[1]);
          eventTime += t0TOF[0] * weight;
          sumOfWeights += weight;
        }

        if (collision.has_foundFT0()) { // T0 measurement is available
          // const auto& ft0 = collision.foundFT0();
          if (collision.t0ACValid()) {
            t0AC[0] = collision.t0AC() * 1000.f;
            t0AC[1] = collision.t0resolution() * 1000.f;
            flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeT0AC;
          }

          weight = 1.f / (t0AC[1] * t0AC[1]);
          eventTime += t0AC[0] * weight;
          sumOfWeights += weight;
        }

        if (sumOfWeights < weightDiamond) { // avoiding sumOfWeights = 0 or worse that diamond
          eventTime = 0;
          sumOfWeights = weightDiamond;
          tableFlags(0);
        } else {
          tableFlags(flags);
        }
        tableEvTime(eventTime / sumOfWeights, std::sqrt(1. / sumOfWeights));
        if (enableTableTOFOnly) {
          tableEvTimeTOFOnly((uint8_t)filterForTOFEventTime(trk), t0TOF[0], t0TOF[1], evTimeTOF.mEventTimeMultiplicity);
        }
      }
    }
  }
  PROCESS_SWITCH(tofEventTime, processFT0, "Process with FT0", false);

  ///
  /// Process function to prepare the event for each track on Run 3 data with only the FT0
  void processOnlyFT0(TrksEvTime const& tracks,
                      aod::FT0s const&,
                      EvTimeCollisionsFT0 const&)
  {
    if (!enableTable) {
      return;
    }

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());
    if (!enableTableTOFOnly) {
      tableEvTimeTOFOnly.reserve(tracks.size());
    }

    for (auto const& t : tracks) { // Loop on collisions
      if (enableTableTOFOnly) {
        tableEvTimeTOFOnly((uint8_t)0, 0.f, 0.f, -1);
      }
      if (!t.has_collision()) { // Track was not assigned, cannot compute event time
        tableFlags(0);
        tableEvTime(0.f, 999.f);
        continue;
      }
      const auto& collision = t.collision_as<EvTimeCollisionsFT0>();

      if (collision.has_foundFT0()) { // T0 measurement is available
        // const auto& ft0 = collision.foundFT0();
        if (collision.t0ACValid()) {
          tableFlags(o2::aod::pidflags::enums::PIDFlags::EvTimeT0AC);
          tableEvTime(collision.t0AC() * 1000.f, collision.t0resolution() * 1000.f);
          continue;
        }
      }
      tableFlags(0);
      tableEvTime(0.f, 999.f);
    }
  }
  PROCESS_SWITCH(tofEventTime, processOnlyFT0, "Process only with FT0", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofSignal>(cfgc)};
  workflow.push_back(adaptAnalysisTask<tofEventTime>(cfgc));
  return workflow;
}
