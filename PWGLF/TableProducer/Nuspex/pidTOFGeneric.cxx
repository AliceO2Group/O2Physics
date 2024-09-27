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
/// \file   pidTOFGeneric.cxx
/// \origin Based on pidTOFBase.cxx
/// \brief  Task to produce event Time obtained from TOF and FT0.
///         In order to redo TOF PID for tracks which are linked to wrong collisions
///

#include <utility>
#include <vector>
#include <string>

// O2 includes
#include "CCDB/BasicCCDBManager.h"
#include "TOFBase/EventTimeMaker.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"

// O2Physics includes
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PID/ParamBase.h"
#include "PID/PIDTOF.h"
#include "PWGLF/DataModel/pidTOFGeneric.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

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
                                                 const float& diamond = 6.0,
                                                 bool isFast = false)
{
  return o2::tof::evTimeMakerFromParam<trackTypeContainer, trackType, trackFilter, response, responseParametersType>(tracks, responseParameters, diamond, isFast);
}

/// Task to produce the event time tables for generic TOF PID
struct pidTOFGeneric {
  // Tables to produce
  Produces<o2::aod::EvTimeTOFFT0> tableEvTime;                 // Table for global event time
  Produces<o2::aod::EvTimeTOFFT0ForTrack> tableEvTimeForTrack; // Table for event time after removing bias from the track
  static constexpr float diamond = 6.0;                        // Collision diamond used in the estimation of the TOF event time
  static constexpr float errDiamond = diamond * 33.356409f;
  static constexpr float weightDiamond = 1.f / (errDiamond * errDiamond);

  bool enableTable = false;
  // Detector response and input parameters
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<bool> inheritFromBaseTask{"inheritFromBaseTask", true, "Flag to iherit all common configurables from the TOF base task"};
  Configurable<bool> fastTOFPID{"fastTOFPID", false, "Flag to enable computeEvTimeFast for evTimeMaker"};
  // CCDB configuration (inherited from TOF signal task)
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
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
    enableTable = isTableRequiredInWorkflow(initContext, "EvTimeTOFFT0") || isTableRequiredInWorkflow(initContext, "EvTimeTOFFT0ForTrack");
    if (!enableTable) {
      LOG(info) << "Table for global Event time is not required, disabling it";
      return;
    }
    LOG(info) << "Table EvTimeTOFFT0 enabled!";

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
  /// Process function to prepare the event time on Run 3 data without the FT0
  using TrksEvTime = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TOFSignal>;
  // Define slice per collision
  Preslice<TrksEvTime> perCollision = aod::track::collisionId;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksEvTime::iterator, pid>;
  using EvTimeCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  void processNoFT0(TrksEvTime& tracks,
                    EvTimeCollisions const& collisions)
  {
    if (!enableTable) {
      return;
    }
    tableEvTime.reserve(collisions.size());
    tableEvTimeForTrack.reserve(tracks.size());

    // for tracks not assigned to a collision
    for (auto track : tracks) {
      if (!track.has_collision()) {
        tableEvTimeForTrack(0.f, 999.f);
      }
    }

    for (auto const& collision : collisions) { // Loop on collisions
      const auto& tracksInCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      if ((sel8TOFEvTime.value == true) && !collision.sel8()) {
        tableEvTime(0.f, 999.f, 0.f, 999.f, 0.f, 999.f);
        for (int i = 0; i < tracksInCollision.size(); i++) {
          tableEvTimeForTrack(0.f, 999.f);
        }
        continue;
      }

      // First make table for event time
      const auto evTimeTOF = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV2, diamond, fastTOFPID);
      int nGoodTracksForTOF = 0; // count for ntrackIndex for removeBias()
      float et = evTimeTOF.mEventTime;
      float erret = evTimeTOF.mEventTimeError;

      if (erret < errDiamond && (maxEvTimeTOF <= 0.f || abs(et) < maxEvTimeTOF)) {
      } else {
        et = 0.f;
        erret = errDiamond;
      }
      tableEvTime(et, erret, et, erret, 0.f, 999.f);
      for (auto const& track : tracksInCollision) {
        evTimeTOF.removeBias<TrksEvTime::iterator, filterForTOFEventTime>(track, nGoodTracksForTOF, et, erret, 2);
        if (erret < errDiamond && (maxEvTimeTOF <= 0.f || abs(et) < maxEvTimeTOF)) {
        } else {
          et = 0.f;
          erret = errDiamond;
        }
        tableEvTimeForTrack(et, erret);
      }
    }
  }
  PROCESS_SWITCH(pidTOFGeneric, processNoFT0, "Process without FT0", false);

  ///
  /// Process function to prepare the event for each track on Run 3 data with the FT0
  using EvTimeCollisionsFT0 = soa::Join<EvTimeCollisions, aod::FT0sCorrected>;
  void processFT0(TrksEvTime& tracks,
                  aod::FT0s const&,
                  EvTimeCollisionsFT0 const& collisions)
  {
    if (!enableTable) {
      return;
    }
    tableEvTime.reserve(collisions.size());
    tableEvTimeForTrack.reserve(tracks.size());
    std::vector<float> tEvTimeForTrack;
    std::vector<float> tEvTimeErrForTrack;
    tEvTimeForTrack.resize(tracks.size());
    tEvTimeErrForTrack.resize(tracks.size());

    // for tracks not assigned to a collision
    for (auto track : tracks) {
      if (!track.has_collision()) {
        tEvTimeForTrack[track.globalIndex()] = 0.f;
        tEvTimeErrForTrack[track.globalIndex()] = 999.f;
      }
    }

    for (auto const& collision : collisions) { // Loop on collisions
      const auto& tracksInCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      if ((sel8TOFEvTime.value == true) && !collision.sel8()) {
        tableEvTime(0.f, 999.f, 0.f, 999.f, 0.f, 999.f);
        for (auto track : tracksInCollision) {
          tEvTimeForTrack[track.globalIndex()] = 0.f;
          tEvTimeErrForTrack[track.globalIndex()] = 999.f;
        }
        continue;
      }

      // Compute the TOF event time
      const auto evTimeTOF = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV2, diamond, fastTOFPID);

      float t0TOF[2] = {static_cast<float_t>(evTimeTOF.mEventTime), static_cast<float_t>(evTimeTOF.mEventTimeError)}; // Value and error of TOF
      float t0AC[2] = {.0f, 999.f};                                                                                   // Value and error of T0A or T0C or T0AC

      float eventTime = 0.f;
      float sumOfWeights = 0.f;
      float weight = 0.f;

      if (t0TOF[1] < errDiamond && (maxEvTimeTOF <= 0 || abs(t0TOF[0]) < maxEvTimeTOF)) {
        weight = 1.f / (t0TOF[1] * t0TOF[1]);
        eventTime += t0TOF[0] * weight;
        sumOfWeights += weight;
      }

      if (collision.has_foundFT0()) { // T0 measurement is available
        // const auto& ft0 = collision.foundFT0();
        if (collision.t0ACValid()) {
          t0AC[0] = collision.t0AC() * 1000.f;
          t0AC[1] = collision.t0resolution() * 1000.f;
        }

        weight = 1.f / (t0AC[1] * t0AC[1]);
        eventTime += t0AC[0] * weight;
        sumOfWeights += weight;
      }

      if (sumOfWeights < weightDiamond) { // avoiding sumOfWeights = 0 or worse that diamond
        eventTime = 0;
        sumOfWeights = weightDiamond;
      }
      tableEvTime(eventTime / sumOfWeights, sqrt(1. / sumOfWeights), t0TOF[0], t0TOF[1], t0AC[0], t0AC[1]);

      int nGoodTracksForTOF = 0; // count for ntrackIndex for removeBias()
      for (auto const& track : tracksInCollision) {
        // Reset the event time
        eventTime = 0.f;
        sumOfWeights = 0.f;
        weight = 0.f;
        // Remove the bias on TOF ev. time
        evTimeTOF.removeBias<TrksEvTime::iterator, filterForTOFEventTime>(track, nGoodTracksForTOF, t0TOF[0], t0TOF[1], 2);
        if (t0TOF[1] < errDiamond && (maxEvTimeTOF <= 0 || abs(t0TOF[0]) < maxEvTimeTOF)) {
          weight = 1.f / (t0TOF[1] * t0TOF[1]);
          eventTime += t0TOF[0] * weight;
          sumOfWeights += weight;
        }

        // Add the contribution from FT0 if it is available, t0AC is already calculated
        if (collision.has_foundFT0()) {
          weight = 1.f / (t0AC[1] * t0AC[1]);
          eventTime += t0AC[0] * weight;
          sumOfWeights += weight;
        }

        if (sumOfWeights < weightDiamond) { // avoiding sumOfWeights = 0 or worse that diamond
          eventTime = 0;
          sumOfWeights = weightDiamond;
        }
        tEvTimeForTrack[track.globalIndex()] = eventTime / sumOfWeights;
        tEvTimeErrForTrack[track.globalIndex()] = sqrt(1. / sumOfWeights);
      }
    }
    for (int i = 0; i < tracks.size(); i++) {
      tableEvTimeForTrack(tEvTimeForTrack[i], tEvTimeErrForTrack[i]);
    }
  }
  PROCESS_SWITCH(pidTOFGeneric, processFT0, "Process with FT0", true);

  ///
  /// Process function to prepare the event time on Run 3 data with only the FT0
  void processOnlyFT0(EvTimeCollisionsFT0 const& collisions,
                      TrksEvTime& tracks,
                      aod::FT0s const&)
  {
    if (!enableTable) {
      return;
    }
    tableEvTime.reserve(collisions.size());
    tableEvTimeForTrack.reserve(tracks.size());

    // for tracks not assigned to a collision
    for (auto track : tracks) {
      if (!track.has_collision()) {
        tableEvTimeForTrack(0.f, 999.f);
      }
    }

    for (auto const& collision : collisions) {
      const auto& tracksInCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      if (collision.has_foundFT0()) { // T0 measurement is available
        // const auto& ft0 = collision.foundFT0();
        if (collision.t0ACValid()) {
          tableEvTime(collision.t0AC() * 1000.f, collision.t0resolution() * 1000.f, 0.f, 999.f, collision.t0AC() * 1000.f, collision.t0resolution() * 1000.f);
          for (int i = 0; i < tracks.size(); i++) {
            tableEvTimeForTrack(collision.t0AC() * 1000.f, collision.t0resolution() * 1000.f);
          }
          return;
        }
      }
      tableEvTime(0.f, 999.f, 0.f, 999.f, 0.f, 999.f);
      for (int i = 0; i < tracksInCollision.size(); i++) {
        tableEvTimeForTrack(0.f, 999.f);
      }
    }
  }
  PROCESS_SWITCH(pidTOFGeneric, processOnlyFT0, "Process only with FT0", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<pidTOFGeneric>(cfgc)};
}
