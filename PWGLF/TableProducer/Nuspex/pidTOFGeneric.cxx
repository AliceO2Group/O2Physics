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
/// \origin Based on pidTOFMerged.cxx
/// \brief  Task to produce event Time obtained from TOF and FT0.
///         In order to redo TOF PID for secondary tracks which are linked to wrong collisions
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
///

#include <string>
#include <utility>
#include <vector>

// O2 includes
#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include "TOFBase/EventTimeMaker.h"

// O2Physics includes
#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PID/PIDTOF.h"
#include "PID/ParamBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

o2::common::core::MetadataHelper metadataInfo;

/// Selection criteria for tracks used for TOF event time
float trackSampleMinMomentum = 0.5f;
float trackSampleMaxMomentum = 2.f;
template <typename trackType>
bool filterForTOFEventTime(const trackType& tr)
{
  return (tr.hasTOF() &&
          tr.p() > trackSampleMinMomentum && tr.p() < trackSampleMaxMomentum &&
          tr.hasITS() &&
          tr.hasTPC() &&
          (tr.trackType() == o2::aod::track::TrackTypeEnum::Track || tr.trackType() == o2::aod::track::TrackTypeEnum::TrackIU));
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
/// Modified based on pidTOFMerge.cxx
struct pidTOFGeneric {
  // Tables to produce
  Produces<o2::aod::EvTimeTOFFT0> tableEvTime;                 // Table for global event time
  Produces<o2::aod::EvTimeTOFFT0ForTrack> tableEvTimeForTrack; // Table for event time after removing bias from the track
  static constexpr bool kRemoveTOFEvTimeBias = true;           // Flag to subtract the Ev. Time bias for low multiplicity events with TOF
  static constexpr float kDiamond = 6.0;                       // Collision diamond used in the estimation of the TOF event time
  static constexpr float kErrDiamond = kDiamond * 33.356409f;
  static constexpr float kWeightDiamond = 1.f / (kErrDiamond * kErrDiamond);

  bool enableTable = false;

  Configurable<bool> fastTOFPID{"fastTOFPID", false, "Flag to enable computeEvTimeFast for evTimeMaker"};
  // Event time configurations
  Configurable<float> minMomentum{"minMomentum", 0.5f, "Minimum momentum to select track sample for TOF event time"};
  Configurable<float> maxMomentum{"maxMomentum", 2.0f, "Maximum momentum to select track sample for TOF event time"};
  Configurable<float> maxEvTimeTOF{"maxEvTimeTOF", 100000.0f, "Maximum value of the TOF event time"};
  Configurable<bool> sel8TOFEvTime{"sel8TOFEvTime", false, "Flag to compute the ev. time only for events that pass the sel8 ev. selection"};
  Configurable<int> mComputeEvTimeWithTOF{"computeEvTimeWithTOF", -1, "Compute ev. time with TOF. -1 (autoset), 0 no, 1 yes"};
  Configurable<int> mComputeEvTimeWithFT0{"computeEvTimeWithFT0", -1, "Compute ev. time with FT0. -1 (autoset), 0 no, 1 yes"};
  Configurable<int> maxNtracksInSet{"maxNtracksInSet", 10, "Size of the set to consider for the TOF ev. time computation"};

  // TOF response and input parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::aod::pidtofgeneric::TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration

  void init(o2::framework::InitContext& initContext)
  {
    mTOFCalibConfig.metadataInfo = metadataInfo;
    mTOFCalibConfig.inheritFromBaseTask(initContext);
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "EvTimeTOFFT0") || isTableRequiredInWorkflow(initContext, "EvTimeTOFFT0ForTrack");
    if (!enableTable) {
      LOG(info) << "Table for global Event time is not required, disabling it";
      // return;  //TODO: uncomment this line
    }
    enableTable = true; // Force enabling the table for now
    LOG(info) << "Table EvTimeTOFFT0 enabled!";

    if (mTOFCalibConfig.autoSetProcessFunctions()) {
      LOG(info) << "Autodetecting process functions";
      if (metadataInfo.isFullyDefined()) {
        if (metadataInfo.isRun3()) {
          doprocessRun3.value = true;
        }
      }
    }

    if (metadataInfo.isFullyDefined()) {
      if (!metadataInfo.isRun3()) {
        LOG(fatal) << "Run3 process not supported in pidTOFGeneric task";
      }
    }

    trackSampleMinMomentum = minMomentum;
    trackSampleMaxMomentum = maxMomentum;
    LOG(info) << "Configuring track sample for TOF ev. time: " << trackSampleMinMomentum << " < p < " << trackSampleMaxMomentum;

    if (sel8TOFEvTime.value == true) {
      LOG(info) << "TOF event time will be computed for collisions that pass the event selection only!";
    }
    mTOFCalibConfig.initSetup(mRespParamsV3, ccdb); // Getting the parametrization parameters

    o2::tof::eventTimeContainer::setMaxNtracksInSet(maxNtracksInSet.value);
    o2::tof::eventTimeContainer::printConfig();
  }

  void process(aod::BCs const&) {}

  ///
  /// Process function to prepare the event for each track on Run 3 data without the FT0
  // Define slice per collision
  using Run3Cols = aod::Collisions;
  using EvTimeCollisions = soa::Join<Run3Cols, aod::EvSels>;
  using EvTimeCollisionsFT0 = soa::Join<EvTimeCollisions, aod::FT0sCorrected>;
  using Run3Trks = o2::soa::Join<aod::TracksIU, aod::TracksExtra>;
  using Run3TrksWtof = soa::Join<Run3Trks, aod::TOFSignal>;
  Preslice<Run3TrksWtof> perCollision = aod::track::collisionId;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<Run3TrksWtof::iterator, pid>;

  void processRun3(Run3TrksWtof const& tracks,
                   aod::FT0s const&,
                   EvTimeCollisionsFT0 const& collisions,
                   aod::BCsWithTimestamps const& bcs)
  {
    if (!enableTable) {
      return;
    }
    LOG(debug) << "Processing Run3 data for TOF event time";

    tableEvTime.reserve(collisions.size());
    tableEvTimeForTrack.reserve(tracks.size());
    mTOFCalibConfig.processSetup(mRespParamsV3, ccdb, bcs.iteratorAt(0)); // Update the calibration parameters

    // Autoset the processing mode for the event time computation
    if (mComputeEvTimeWithTOF == -1 || mComputeEvTimeWithFT0 == -1) {
      switch (mTOFCalibConfig.collisionSystem()) {
        case CollisionSystemType::kCollSyspp: // pp
          mComputeEvTimeWithTOF.value = ((mComputeEvTimeWithTOF == -1) ? 0 : mComputeEvTimeWithTOF.value);
          mComputeEvTimeWithFT0.value = ((mComputeEvTimeWithFT0 == -1) ? 1 : mComputeEvTimeWithFT0.value);
          break;
        case CollisionSystemType::kCollSysPbPb: // PbPb
          mComputeEvTimeWithTOF.value = ((mComputeEvTimeWithTOF == -1) ? 1 : mComputeEvTimeWithTOF.value);
          mComputeEvTimeWithFT0.value = ((mComputeEvTimeWithFT0 == -1) ? 0 : mComputeEvTimeWithFT0.value);
          break;
        default:
          LOG(fatal) << "Collision system " << mTOFCalibConfig.collisionSystem() << " " << CollisionSystemType::getCollisionSystemName(mTOFCalibConfig.collisionSystem()) << " not supported for TOF event time computation";
          break;
      }
    }
    LOG(debug) << "Running on " << CollisionSystemType::getCollisionSystemName(mTOFCalibConfig.collisionSystem()) << " mComputeEvTimeWithTOF " << mComputeEvTimeWithTOF.value << " mComputeEvTimeWithFT0 " << mComputeEvTimeWithFT0.value;

    if (mComputeEvTimeWithTOF == 1 && mComputeEvTimeWithFT0 == 1) {
      int lastCollisionId = -1;                                                                                       // Last collision ID analysed
      for (auto const& t : tracks) {                                                                                  // Loop on collisions
        if (!t.has_collision() || ((sel8TOFEvTime.value == true) && !t.collision_as<EvTimeCollisionsFT0>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
          tableEvTimeForTrack(0.f, 999.f);
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
        const auto evTimeMakerTOF = evTimeMakerForTracks<Run3TrksWtof::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV3, kDiamond, fastTOFPID);

        float t0AC[2] = {.0f, 999.f};                                                                                             // Value and error of T0A or T0C or T0AC
        float t0TOF[2] = {static_cast<float_t>(evTimeMakerTOF.mEventTime), static_cast<float_t>(evTimeMakerTOF.mEventTimeError)}; // Value and error of TOF

        int nGoodTracksForTOF = 0;
        float eventTime = 0.f;
        float sumOfWeights = 0.f;
        float weight = 0.f;

        if (t0TOF[1] < kErrDiamond && (maxEvTimeTOF <= 0 || std::abs(t0TOF[0]) < maxEvTimeTOF)) {
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

        tableEvTime(eventTime / sumOfWeights, std::sqrt(1. / sumOfWeights), t0TOF[0], t0TOF[1], t0AC[0], t0AC[1]);

        for (auto const& trk : tracksInCollision) { // Loop on Tracks
          // Reset the event time
          eventTime = 0.f;
          sumOfWeights = 0.f;
          weight = 0.f;
          // Remove the bias on TOF ev. time
          if constexpr (kRemoveTOFEvTimeBias) {
            evTimeMakerTOF.removeBias<Run3TrksWtof::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, t0TOF[0], t0TOF[1], 2);
          }
          if (t0TOF[1] < kErrDiamond && (maxEvTimeTOF <= 0 || std::abs(t0TOF[0]) < maxEvTimeTOF)) {
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

          if (sumOfWeights < kWeightDiamond) { // avoiding sumOfWeights = 0 or worse that kDiamond
            eventTime = 0;
            sumOfWeights = kWeightDiamond;
          }
          tableEvTimeForTrack(eventTime / sumOfWeights, std::sqrt(1. / sumOfWeights));
        }
      }
    } else if (mComputeEvTimeWithTOF == 1 && mComputeEvTimeWithFT0 == 0) {
      int lastCollisionId = -1;                                                                                    // Last collision ID analysed
      for (auto const& t : tracks) {                                                                               // Loop on collisions
        if (!t.has_collision() || ((sel8TOFEvTime.value == true) && !t.collision_as<EvTimeCollisions>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
          tableEvTimeForTrack(0.f, 999.f);
          continue;
        }
        if (t.collisionId() == lastCollisionId) { // Event time from this collision is already in the table
          continue;
        }
        /// Create new table for the tracks in a collision
        lastCollisionId = t.collisionId(); /// Cache last collision ID

        const auto& tracksInCollision = tracks.sliceBy(perCollision, lastCollisionId);

        // First make table for event time
        const auto evTimeMakerTOF = evTimeMakerForTracks<Run3TrksWtof::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV3, kDiamond, fastTOFPID);
        int nGoodTracksForTOF = 0;
        float et = evTimeMakerTOF.mEventTime;
        float erret = evTimeMakerTOF.mEventTimeError;

        if (erret < kErrDiamond && (maxEvTimeTOF <= 0.f || std::abs(et) < maxEvTimeTOF)) {
        } else {
          et = 0.f;
          erret = kErrDiamond;
        }
        tableEvTime(et, erret, et, erret, 0.f, 999.f);

        for (auto const& trk : tracksInCollision) { // Loop on Tracks
          if constexpr (kRemoveTOFEvTimeBias) {
            evTimeMakerTOF.removeBias<Run3TrksWtof::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, et, erret, 2);
          }
          if (erret < kErrDiamond && (maxEvTimeTOF <= 0.f || std::abs(et) < maxEvTimeTOF)) {
          } else {
            et = 0.f;
            erret = kErrDiamond;
          }
          tableEvTimeForTrack(et, erret);
        }
      }
    } else if (mComputeEvTimeWithTOF == 0 && mComputeEvTimeWithFT0 == 1) {
      for (const auto& track : tracks) {
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
            for (int i = 0; i < tracksInCollision.size(); i++) {
              tableEvTimeForTrack(collision.t0AC() * 1000.f, collision.t0resolution() * 1000.f);
            }
            continue;
          }
        }
        tableEvTime(0.f, 999.f, 0.f, 999.f, 0.f, 999.f);
        for (int i = 0; i < tracksInCollision.size(); i++) {
          tableEvTimeForTrack(0.f, 999.f);
        }
      }
    } else {
      LOG(fatal) << "Invalid configuration for TOF event time computation";
    }
  }
  PROCESS_SWITCH(pidTOFGeneric, processRun3, "Process the Run3 data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{
    adaptAnalysisTask<pidTOFGeneric>(cfgc)};
}
