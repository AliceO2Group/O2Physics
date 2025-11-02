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
/// \file pidTOFMerge.cxx
///
/// \brief Task to produce PID tables for TOF split for each particle.
///        Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
///

#include "pidTOFBase.h"

#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/Core/PID/PIDTOFParamService.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseTOF.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <DataFormatsTOF/ParameterContainers.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
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

#include <TGraph.h>
#include <TH2.h>
#include <TString.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

// Input data types
using Run3Trks = o2::soa::Join<aod::TracksIU, aod::TracksExtra>;
using Run3TrksWtof = soa::Join<Run3Trks, aod::TOFSignal>;
using Run3TrksWtofWevTime = soa::Join<Run3TrksWtof, aod::TOFEvTime, aod::pidEvTimeFlags>;

using EvTimeCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using EvTimeCollisionsFT0 = soa::Join<EvTimeCollisions, aod::FT0sCorrected>;

using Run2Trks = o2::soa::Join<aod::Tracks, aod::TracksExtra>;
using Run2TrksWtofWevTime = soa::Join<Run2Trks, aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags>;

/// Selection criteria for tracks used for TOF event time
bool isTrackGoodMatchForTOFPID(const Run3Trks::iterator& tr)
{
  if (!tr.hasTOF()) {
    return false;
  }
  return true;
}

/// Task to produce the TOF signal from the trackTime information
struct tofSignal {
  // Detector response and input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::pid::tof::TOFResponse> tofResponse;
  // Tables to produce
  o2::framework::Produces<o2::aod::TOFSignal> table;
  o2::framework::Produces<o2::aod::pidTOFFlags> tableFlags;
  // Running flags
  bool enableTableTOFSignal = false;   // Flag to check if the TOF signal table is requested or not
  bool enableTablepidTOFFlags = false; // Flag to check if the TOF signal flags table is requested or not
  // Output histograms
  Configurable<bool> enableQaHistograms{"enableQaHistograms", false, "Flag to enable the QA histograms"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  struct : ConfigurableGroup {
    Configurable<std::string> cfgUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> cfgPathGrpLhcIf{"ccdb-path-grplhcif", "GLO/Config/GRPLHCIF", "Path on the CCDB for the GRPLHCIF object"};
    Configurable<int64_t> cfgTimestamp{"ccdb-timestamp", -1, "timestamp of the object"};
    Configurable<std::string> cfgTimeShiftCCDBPathPos{"timeShiftCCDBPathPos", "", "Path of the TOF time shift vs eta for pos. tracks. If empty none is taken"};
    Configurable<std::string> cfgTimeShiftCCDBPathNeg{"timeShiftCCDBPathNeg", "", "Path of the TOF time shift vs eta for neg. tracks. If empty none is taken"};
    Configurable<std::string> cfgTimeShiftCCDBPathPosMC{"timeShiftCCDBPathPosMC", "", "Path of the TOF time shift for MC vs eta for pos. tracks. If empty none is taken"};
    Configurable<std::string> cfgTimeShiftCCDBPathNegMC{"timeShiftCCDBPathNegMC", "", "Path of the TOF time shift for MC vs eta for neg. tracks. If empty none is taken"};
    Configurable<std::string> cfgParamFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
    Configurable<std::string> cfgParametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
    Configurable<std::string> cfgReconstructionPass{"reconstructionPass", "", {"Apass to use when fetching the calibration tables. Empty (default) does not check for any pass. Use `metadata` to fetch it from the AO2D metadata. Otherwise it will override the metadata."}};
    Configurable<std::string> cfgReconstructionPassDefault{"reconstructionPassDefault", "unanchored", {"Default pass to get if the standard one is not found"}};
    Configurable<bool> cfgFatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};
    Configurable<bool> cfgEnableTimeDependentResponse{"enableTimeDependentResponse", false, "Flag to use the collision timestamp to fetch the PID Response"};
    Configurable<int> cfgCollisionSystem{"collisionSystem", -1, "Collision system: -1 (autoset), 0 (pp), 1 (PbPb), 2 (XeXe), 3 (pPb)"};
    Configurable<bool> cfgAutoSetProcessFunctions{"autoSetProcessFunctions", true, "Flag to autodetect the process functions to use"};
  } cfg; // Configurables (only defined here and inherited from other tasks)

  void init(o2::framework::InitContext& initContext)
  {
    LOG(debug) << "Initializing the tofSignal task";
    tofResponse->initSetup(ccdb, initContext);
    // Checking that the table is requested in the workflow and enabling it
    enableTableTOFSignal = isTableRequiredInWorkflow(initContext, "TOFSignal");
    if (enableTableTOFSignal) {
      LOG(info) << "Table TOFSignal enabled!";
    }
    enableTablepidTOFFlags = isTableRequiredInWorkflow(initContext, "pidTOFFlags");
    if (enableTablepidTOFFlags) {
      LOG(info) << "Table pidTOFFlags enabled!";
    }

    // If the table is not requested, disable the task. Uless a process function is enabled from the workflow configuration
    if (!enableTableTOFSignal && !enableTablepidTOFFlags && !doprocessRun2 && !doprocessRun3) {
      LOG(info) << "No table or process is enabled. Disabling task";
      return;
    }
    if (tofResponse->cfgAutoSetProcessFunctions()) {
      LOG(info) << "Autodetecting process functions";
      if (tofResponse->metadataInfo.isFullyDefined() && !doprocessRun2 && !doprocessRun3) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
        if (tofResponse->metadataInfo.isRun3()) {
          doprocessRun3.value = true;
        } else {
          doprocessRun2.value = false;
        }
      }
    }

    // Last checks on the process functions
    if (doprocessRun2 && doprocessRun3) {
      LOG(fatal) << "Both processRun2 and processRun3 are enabled. Pick one of the two";
    }
    if (!doprocessRun2 && !doprocessRun3) {
      LOG(fatal) << "Neither processRun2 nor processRun3 are enabled. Pick one of the two";
    }
    if (!enableQaHistograms) {
      return;
    }
    histos.add("tofSignal", "tofSignal", kTH1D, {{1000, -1000, 1000000, "tofSignal (ps)"}});
    if (enableTablepidTOFFlags) {
      histos.add("goodForPIDFlags", "goodForPIDFlags", kTH1D, {{3, 0, 3, "flags"}});
    }
  }

  /// Dummy process function for BCs, needed in case both Run2 and Run3 process functions are disabled
  void process(aod::BCs const&) {}

  void processRun3(Run3Trks const& tracks)
  {
    if (!enableTableTOFSignal) {
      return;
    }
    table.reserve(tracks.size());
    if (enableTablepidTOFFlags) {
      tableFlags.reserve(tracks.size());
    }
    for (const auto& trk : tracks) {
      const float& sig = o2::pid::tof::TOFSignal<Run3Trks::iterator>::GetTOFSignal(trk);
      if (enableQaHistograms) {
        histos.fill(HIST("tofSignal"), sig);
      }
      table(sig);
      if (!enableTablepidTOFFlags) {
        continue;
      }
      const auto& b = isTrackGoodMatchForTOFPID(trk);
      if (enableQaHistograms) {
        histos.fill(HIST("goodForPIDFlags"), sig);
      }
      tableFlags(b);
    }
  }
  PROCESS_SWITCH(tofSignal, processRun3, "Process Run3 data i.e. input is TrackIU. Set to false to autodetect from metadata.", false);

  void processRun2(Run2Trks const& tracks)
  {
    if (!enableTableTOFSignal) {
      return;
    }
    table.reserve(tracks.size());
    if (enableTablepidTOFFlags) {
      tableFlags.reserve(tracks.size());
    }
    for (const auto& trk : tracks) {
      table(o2::pid::tof::TOFSignal<Run2Trks::iterator>::GetTOFSignal(trk));
      if (!enableTablepidTOFFlags) {
        continue;
      }
      tableFlags(true);
    }
  }
  PROCESS_SWITCH(tofSignal, processRun2, "Process Run2 data i.e. input is Tracks. Set to false to autodetect from metadata.", false);
};

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
                                                 const float& diamond = 6.0)
{
  return o2::tof::evTimeMakerFromParam<trackTypeContainer, trackType, trackFilter, response, responseParametersType>(tracks, responseParameters, diamond);
}

// Part 2 event time definition

/// Task to produce the TOF event time table
struct tofEventTime {
  // Detector response and input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::pid::tof::TOFResponse> tofResponse;
  // Tables to produce
  Produces<o2::aod::TOFEvTime> tableEvTime;
  Produces<o2::aod::EvTimeTOFOnly> tableEvTimeTOFOnly;
  Produces<o2::aod::pidEvTimeFlags> tableFlags;
  static constexpr bool kRemoveTOFEvTimeBias = true; // Flag to subtract the Ev. Time bias for low multiplicity events with TOF
  static constexpr float kDiamond = 6.0;             // Collision diamond used in the estimation of the TOF event time
  static constexpr float kErrDiamond = kDiamond * 33.356409f;
  static constexpr float kWeightDiamond = 1.f / (kErrDiamond * kErrDiamond);

  bool enableTableTOFEvTime = false;
  bool enableTableEvTimeTOFOnly = false;

  // Event time configurations
  Configurable<float> minMomentum{"minMomentum", 0.5f, "Minimum momentum to select track sample for TOF event time"};
  Configurable<float> maxMomentum{"maxMomentum", 2.0f, "Maximum momentum to select track sample for TOF event time"};
  Configurable<float> maxEvTimeTOF{"maxEvTimeTOF", 100000.0f, "Maximum value of the TOF event time"};
  Configurable<bool> sel8TOFEvTime{"sel8TOFEvTime", false, "Flag to compute the ev. time only for events that pass the sel8 ev. selection"};
  Configurable<int> mComputeEvTimeWithTOF{"computeEvTimeWithTOF", -1, "Compute ev. time with TOF. -1 (autoset), 0 no, 1 yes"};
  Configurable<int> mComputeEvTimeWithFT0{"computeEvTimeWithFT0", -1, "Compute ev. time with FT0. -1 (autoset), 0 no, 1 yes"};
  Configurable<int> maxNtracksInSet{"maxNtracksInSet", 10, "Size of the set to consider for the TOF ev. time computation"};

  void init(o2::framework::InitContext& initContext)
  {
    LOG(debug) << "Initializing the tofEventTime task";
    tofResponse->initSetup(ccdb, initContext);
    // Checking that the table is requested in the workflow and enabling it
    enableTableTOFEvTime = isTableRequiredInWorkflow(initContext, "TOFEvTime");

    if (!enableTableTOFEvTime) {
      LOG(info) << "Table for TOF Event time (TOFEvTime) is not required, disabling it";
    }
    LOG(info) << "Table TOFEvTime enabled!";

    enableTableEvTimeTOFOnly = isTableRequiredInWorkflow(initContext, "EvTimeTOFOnly");
    if (enableTableEvTimeTOFOnly) {
      LOG(info) << "Table EvTimeTOFOnly enabled!";
    }

    if (!enableTableTOFEvTime && !enableTableEvTimeTOFOnly) {
      LOG(info) << "No table is enabled. Disabling task";
      return;
    }

    if (tofResponse->cfgAutoSetProcessFunctions()) {
      LOG(info) << "Autodetecting process functions";
      if (tofResponse->metadataInfo.isFullyDefined()) {
        if (tofResponse->metadataInfo.isRun3()) {
          doprocessRun3.value = true;
        } else {
          doprocessRun2.value = true;
        }
      }
    }

    if (tofResponse->metadataInfo.isFullyDefined()) {
      if (tofResponse->metadataInfo.isRun3() && doprocessRun2) {
        LOG(fatal) << "Run2 process function is enabled but the metadata says it is Run3";
      }
      if (!tofResponse->metadataInfo.isRun3() && doprocessRun3) {
        LOG(fatal) << "Run3 process function is enabled but the metadata says it is Run2";
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
    if (doprocessRun3 == true) {
      LOGF(info, "Enabling process function: processRun3");
      nEnabled++;
    }
    if (nEnabled > 1) {
      LOGF(fatal, "Cannot enable more process functions at the same time. Please choose one.");
    }

    if (sel8TOFEvTime.value == true) {
      LOG(info) << "TOF event time will be computed for collisions that pass the event selection only!";
    }
    o2::tof::eventTimeContainer::setMaxNtracksInSet(maxNtracksInSet.value);
    o2::tof::eventTimeContainer::printConfig();
  }

  void process(aod::BCs const&) {}

  ///
  /// Process function to prepare the event for each track on Run 2 data
  void processRun2(aod::Tracks const& tracks,
                   aod::Collisions const&,
                   aod::BCsWithTimestamps const& bcs)
  {
    if (!enableTableTOFEvTime) {
      return;
    }
    tofResponse->processSetup(bcs.iteratorAt(0)); // Update the response parameters

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
  PROCESS_SWITCH(tofEventTime, processRun2, "Process with Run2 data", true);

  ///
  /// Process function to prepare the event for each track on Run 3 data without the FT0
  // Define slice per collision
  Preslice<Run3TrksWtof> perCollision = aod::track::collisionId;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<Run3TrksWtof::iterator, pid>;
  void processRun3(Run3TrksWtof const& tracks,
                   aod::FT0s const&,
                   EvTimeCollisionsFT0 const&,
                   aod::BCsWithTimestamps const& bcs)
  {
    if (!enableTableTOFEvTime) {
      return;
    }
    LOG(debug) << "Processing Run3 data for TOF event time";

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());
    if (enableTableEvTimeTOFOnly) {
      tableEvTimeTOFOnly.reserve(tracks.size());
    }
    tofResponse->processSetup(bcs.iteratorAt(0)); // Update the response parameters

    // Autoset the processing mode for the event time computation
    if (mComputeEvTimeWithTOF == -1 || mComputeEvTimeWithFT0 == -1) {
      switch (tofResponse->cfgCollisionType()) {
        case CollisionSystemType::kCollSyspp: // pp
          mComputeEvTimeWithTOF.value = ((mComputeEvTimeWithTOF == -1) ? 0 : mComputeEvTimeWithTOF.value);
          mComputeEvTimeWithFT0.value = ((mComputeEvTimeWithFT0 == -1) ? 1 : mComputeEvTimeWithFT0.value);
          break;
        case CollisionSystemType::kCollSysPbPb: // PbPb
          mComputeEvTimeWithTOF.value = ((mComputeEvTimeWithTOF == -1) ? 1 : mComputeEvTimeWithTOF.value);
          mComputeEvTimeWithFT0.value = ((mComputeEvTimeWithFT0 == -1) ? 0 : mComputeEvTimeWithFT0.value);
          break;
        default:
          LOG(fatal) << "Collision system " << tofResponse->cfgCollisionType() << " " << CollisionSystemType::getCollisionSystemName(tofResponse->cfgCollisionType()) << " not supported for TOF event time computation";
          break;
      }
    }
    LOG(debug) << "Running on " << CollisionSystemType::getCollisionSystemName(tofResponse->cfgCollisionType()) << " mComputeEvTimeWithTOF " << mComputeEvTimeWithTOF.value << " mComputeEvTimeWithFT0 " << mComputeEvTimeWithFT0.value;

    if (mComputeEvTimeWithTOF == 1 && mComputeEvTimeWithFT0 == 1) {
      int lastCollisionId = -1;                                                                                       // Last collision ID analysed
      for (auto const& t : tracks) {                                                                                  // Loop on collisions
        if (!t.has_collision() || ((sel8TOFEvTime.value == true) && !t.collision_as<EvTimeCollisionsFT0>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
          tableFlags(0);
          tableEvTime(0.f, 999.f);
          if (enableTableEvTimeTOFOnly) {
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
        const auto evTimeMakerTOF = evTimeMakerForTracks<Run3TrksWtof::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, tofResponse->parameters, kDiamond);
        float t0AC[2] = {.0f, 999.f};                                                                                             // Value and error of T0A or T0C or T0AC
        float t0TOF[2] = {static_cast<float_t>(evTimeMakerTOF.mEventTime), static_cast<float_t>(evTimeMakerTOF.mEventTimeError)}; // Value and error of TOF

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
          if constexpr (kRemoveTOFEvTimeBias) {
            evTimeMakerTOF.removeBias<Run3TrksWtof::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, t0TOF[0], t0TOF[1], 2);
          }
          if (t0TOF[1] < kErrDiamond && (maxEvTimeTOF <= 0 || std::abs(t0TOF[0]) < maxEvTimeTOF)) {
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

          if (sumOfWeights < kWeightDiamond) { // avoiding sumOfWeights = 0 or worse that kDiamond
            eventTime = 0;
            sumOfWeights = kWeightDiamond;
            tableFlags(0);
          } else {
            tableFlags(flags);
          }
          tableEvTime(eventTime / sumOfWeights, std::sqrt(1. / sumOfWeights));
          if (enableTableEvTimeTOFOnly) {
            tableEvTimeTOFOnly((uint8_t)filterForTOFEventTime(trk), t0TOF[0], t0TOF[1], evTimeMakerTOF.mEventTimeMultiplicity);
          }
        }
      }
    } else if (mComputeEvTimeWithTOF == 1 && mComputeEvTimeWithFT0 == 0) {
      int lastCollisionId = -1;                                                                                    // Last collision ID analysed
      for (auto const& t : tracks) {                                                                               // Loop on collisions
        if (!t.has_collision() || ((sel8TOFEvTime.value == true) && !t.collision_as<EvTimeCollisions>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
          tableFlags(0);
          tableEvTime(0.f, 999.f);
          if (enableTableEvTimeTOFOnly) {
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
        const auto evTimeMakerTOF = evTimeMakerForTracks<Run3TrksWtof::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, tofResponse->parameters, kDiamond);
        int nGoodTracksForTOF = 0;
        float et = evTimeMakerTOF.mEventTime;
        float erret = evTimeMakerTOF.mEventTimeError;

        for (auto const& trk : tracksInCollision) { // Loop on Tracks
          if constexpr (kRemoveTOFEvTimeBias) {
            evTimeMakerTOF.removeBias<Run3TrksWtof::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, et, erret, 2);
          }
          uint8_t flags = 0;
          if (erret < kErrDiamond && (maxEvTimeTOF <= 0.f || std::abs(et) < maxEvTimeTOF)) {
            flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeTOF;
          } else {
            et = 0.f;
            erret = kErrDiamond;
          }
          tableFlags(flags);
          tableEvTime(et, erret);
          if (enableTableEvTimeTOFOnly) {
            tableEvTimeTOFOnly((uint8_t)filterForTOFEventTime(trk), et, erret, evTimeMakerTOF.mEventTimeMultiplicity);
          }
        }
      }
    } else if (mComputeEvTimeWithTOF == 0 && mComputeEvTimeWithFT0 == 1) {
      for (auto const& t : tracks) { // Loop on collisions
        if (enableTableEvTimeTOFOnly) {
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
    } else {
      LOG(fatal) << "Invalid configuration for TOF event time computation";
    }
  }
  PROCESS_SWITCH(tofEventTime, processRun3, "Process the Run3 data", true);
};

// Part 3 Nsigma computation

static constexpr int kParEnabledN = 2;
static constexpr int kIdxEl = 0;
static constexpr int kIdxMu = 1;
static constexpr int kIdxPi = 2;
static constexpr int kIdxKa = 3;
static constexpr int kIdxPr = 4;
static constexpr int kIdxDe = 5;
static constexpr int kIdxTr = 6;
static constexpr int kIdxHe = 7;
static constexpr int kIdxAl = 8;

static const std::vector<std::string> kParEnabledNames{"Enable", "EnableFull"};
static constexpr int kDefaultParEnabled[nSpecies][kParEnabledN]{{-1, -1},
                                                                {-1, -1},
                                                                {-1, -1},
                                                                {-1, -1},
                                                                {-1, -1},
                                                                {-1, -1},
                                                                {-1, -1},
                                                                {-1, -1},
                                                                {-1, -1}};

/// Task to produce the response table
struct tofPidMerge {
  // Detector response and input parameters
  Service<o2::pid::tof::TOFResponse> tofResponse;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Tables to produce
  Produces<o2::aod::pidTOFEl> tablePIDEl;
  Produces<o2::aod::pidTOFMu> tablePIDMu;
  Produces<o2::aod::pidTOFPi> tablePIDPi;
  Produces<o2::aod::pidTOFKa> tablePIDKa;
  Produces<o2::aod::pidTOFPr> tablePIDPr;
  Produces<o2::aod::pidTOFDe> tablePIDDe;
  Produces<o2::aod::pidTOFTr> tablePIDTr;
  Produces<o2::aod::pidTOFHe> tablePIDHe;
  Produces<o2::aod::pidTOFAl> tablePIDAl;

  // Tables to produce (full)
  Produces<o2::aod::pidTOFFullEl> tablePIDFullEl;
  Produces<o2::aod::pidTOFFullMu> tablePIDFullMu;
  Produces<o2::aod::pidTOFFullPi> tablePIDFullPi;
  Produces<o2::aod::pidTOFFullKa> tablePIDFullKa;
  Produces<o2::aod::pidTOFFullPr> tablePIDFullPr;
  Produces<o2::aod::pidTOFFullDe> tablePIDFullDe;
  Produces<o2::aod::pidTOFFullTr> tablePIDFullTr;
  Produces<o2::aod::pidTOFFullHe> tablePIDFullHe;
  Produces<o2::aod::pidTOFFullAl> tablePIDFullAl;

  // Beta tables
  Produces<aod::pidTOFbeta> tablePIDBeta;
  Produces<aod::pidTOFmass> tablePIDTOFMass;
  bool enableTableBeta = false;
  bool enableTableMass = false;

  Configurable<bool> enableQaHistograms{"enableQaHistograms", false, "Flag to enable the QA histograms"};
  Configurable<bool> enableTOFParamsForBetaMass{"enableTOFParamsForBetaMass", false, "Flag to use TOF parameters for TOF Beta and Mass"};

  // Configuration flags to include and exclude particle hypotheses
  Configurable<LabeledArray<int>> enableParticle{"enableParticle",
                                                 {kDefaultParEnabled[0], nSpecies, kParEnabledN, particleNames, kParEnabledNames},
                                                 "Produce PID information for the various mass hypotheses. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};

  // Histograms for QA
  std::array<std::shared_ptr<TH2>, nSpecies> hnsigma;
  std::array<std::shared_ptr<TH2>, nSpecies> hnsigmaFull;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Running variables
  std::vector<int> mEnabledParticles;     // Vector of enabled PID hypotheses to loop on when making tables
  std::vector<int> mEnabledParticlesFull; // Vector of enabled PID hypotheses to loop on when making full tables
  void init(o2::framework::InitContext& initContext)
  {
    LOG(debug) << "Initializing the TOF PID Merge task";
    tofResponse->initSetup(ccdb, initContext);
    // Checking the tables are requested in the workflow and enabling them
    for (int i = 0; i < nSpecies; i++) {
      // First checking tiny
      int f = enableParticle->get(particleNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(initContext, "pidTOF" + particleNames[i], f);
      if (f == 1) {
        mEnabledParticles.push_back(i);
      }

      // Then checking full tables
      f = enableParticle->get(particleNames[i].c_str(), "EnableFull");
      enableFlagIfTableRequired(initContext, "pidTOFFull" + particleNames[i], f);
      if (f == 1) {
        mEnabledParticlesFull.push_back(i);
      }
    }
    if (mEnabledParticlesFull.size() == 0 && mEnabledParticles.size() == 0) {
      LOG(info) << "No PID tables are required, disabling the task";
      doprocessRun3.value = false;
      doprocessRun2.value = false;
    } else {
      if (tofResponse->cfgAutoSetProcessFunctions()) {
        LOG(info) << "Autodetecting process functions for mass and beta";
        if (tofResponse->metadataInfo.isFullyDefined()) {
          if (tofResponse->metadataInfo.isRun3()) {
            doprocessRun3.value = true;
            doprocessRun2.value = false;
          } else {
            doprocessRun2.value = true;
            doprocessRun3.value = false;
          }
        }
      }
      if (doprocessRun2 && doprocessRun3) {
        LOG(fatal) << "Both processRun2 and processRun3 are enabled. Pick one of the two";
      }
      if (!doprocessRun2 && !doprocessRun3) {
        LOG(fatal) << "Neither processRun2 nor processRun3 are enabled. Pick one of the two";
      }
    }

    // Printing enabled tables and enabling QA histograms if needed
    LOG(info) << "++ Enabled tables:";
    const AxisSpec pAxis{100, 0, 5, "#it{p} (GeV/#it{c})"};
    const AxisSpec nSigmaAxis{100, -10, 10, "N_{#sigma}^{TOF}"};
    for (const int& i : mEnabledParticles) {
      LOG(info) << "++  pidTOF" << particleNames[i] << " is enabled";
      if (!enableQaHistograms) {
        continue;
      }
      hnsigma[i] = histos.add<TH2>(Form("nsigma/%s", particleNames[i].c_str()), Form("N_{#sigma}^{TOF}(%s)", particleNames[i].c_str()), kTH2F, {pAxis, nSigmaAxis});
    }
    for (const int& i : mEnabledParticlesFull) {
      LOG(info) << "++  pidTOFFull" << particleNames[i] << " is enabled";
      if (!enableQaHistograms) {
        continue;
      }
      hnsigmaFull[i] = histos.add<TH2>(Form("nsigmaFull/%s", particleNames[i].c_str()), Form("N_{#sigma}^{TOF}(%s)", particleNames[i].c_str()), kTH2F, {pAxis, nSigmaAxis});
    }

    // Checking the TOF mass and TOF beta tables
    enableTableBeta = isTableRequiredInWorkflow(initContext, "pidTOFbeta");
    enableTableMass = isTableRequiredInWorkflow(initContext, "pidTOFmass");

    if (!enableTableBeta && !enableTableMass) {
      LOG(info) << "No table for TOF mass and beta is required. Disabling beta and mass tables";
      doprocessRun2BetaM.value = false;
      doprocessRun3BetaM.value = false;
    } else {
      if (tofResponse->cfgAutoSetProcessFunctions()) {
        LOG(info) << "Autodetecting process functions for mass and beta";
        if (tofResponse->metadataInfo.isFullyDefined()) {
          if (tofResponse->metadataInfo.isRun3()) {
            doprocessRun3BetaM.value = true;
            doprocessRun2BetaM.value = false;
          } else {
            doprocessRun2BetaM.value = true;
            doprocessRun3BetaM.value = false;
          }
        } else {
          tofResponse->metadataInfo.print();
          LOG(warning) << "Metadata is not defined, cannot autodetect process functions for mass and beta";
        }
      } else {
        LOG(info) << "Process functions for mass and beta are set manually";
      }
      if (doprocessRun2BetaM && doprocessRun3BetaM) {
        LOG(fatal) << "Both processRun2BetaM and processRun3BetaM are enabled. Pick one of the two";
      }
      if (!doprocessRun2BetaM && !doprocessRun3BetaM) {
        LOG(fatal) << "Neither processRun2BetaM nor processRun3BetaM are enabled. Pick one of the two";
      }
    }
  }

  // Reserves an empty table for the given particle ID with size of the given track table
  void reserveTable(const int id, const int64_t& size, const bool fullTable = false)
  {
    switch (id) {
      case kIdxEl: {
        if (fullTable) {
          tablePIDFullEl.reserve(size);
        } else {
          tablePIDEl.reserve(size);
        }
        break;
      }
      case kIdxMu: {
        if (fullTable) {
          tablePIDFullMu.reserve(size);
        } else {
          tablePIDMu.reserve(size);
        }
        break;
      }
      case kIdxPi: {
        if (fullTable) {
          tablePIDFullPi.reserve(size);
        } else {
          tablePIDPi.reserve(size);
        }
        break;
      }
      case kIdxKa: {
        if (fullTable) {
          tablePIDFullKa.reserve(size);
        } else {
          tablePIDKa.reserve(size);
        }
        break;
      }
      case kIdxPr: {
        if (fullTable) {
          tablePIDFullPr.reserve(size);
        } else {
          tablePIDPr.reserve(size);
        }
        break;
      }
      case kIdxDe: {
        if (fullTable) {
          tablePIDFullDe.reserve(size);
        } else {
          tablePIDDe.reserve(size);
        }
        break;
      }
      case kIdxTr: {
        if (fullTable) {
          tablePIDFullTr.reserve(size);
        } else {
          tablePIDTr.reserve(size);
        }
        break;
      }
      case kIdxHe: {
        if (fullTable) {
          tablePIDFullHe.reserve(size);
        } else {
          tablePIDHe.reserve(size);
        }
        break;
      }
      case kIdxAl: {
        if (fullTable) {
          tablePIDFullAl.reserve(size);
        } else {
          tablePIDAl.reserve(size);
        }
        break;
      }
      default:
        LOG(fatal) << "Wrong particle ID in reserveTable() for " << (fullTable ? "full" : "tiny") << " tables";
        break;
    }
  }

  // Makes the table empty for the given particle ID, filling it with dummy values
  void makeTableEmpty(const int id, bool fullTable = false)
  {
    switch (id) {
      case kIdxEl:
        if (fullTable) {
          tablePIDFullEl(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDEl);
        }
        break;
      case kIdxMu:
        if (fullTable) {
          tablePIDFullMu(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDMu);
        }
        break;
      case kIdxPi:
        if (fullTable) {
          tablePIDFullPi(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDPi);
        }
        break;
      case kIdxKa:
        if (fullTable) {
          tablePIDFullKa(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDKa);
        }
        break;
      case kIdxPr:
        if (fullTable) {
          tablePIDFullPr(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDPr);
        }
        break;
      case kIdxDe:
        if (fullTable) {
          tablePIDFullDe(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDDe);
        }
        break;
      case kIdxTr:
        if (fullTable) {
          tablePIDFullTr(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDTr);
        }
        break;
      case kIdxHe:
        if (fullTable) {
          tablePIDFullHe(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDHe);
        }
        break;
      case kIdxAl:
        if (fullTable) {
          tablePIDFullAl(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDAl);
        }
        break;
      default:
        LOG(fatal) << "Wrong particle ID in makeTableEmpty() for " << (fullTable ? "full" : "tiny") << " tables";
        break;
    }
  }

  void process(aod::BCs const&) {}

  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Run3TrksWtofWevTime::iterator, pid>;
  void processRun3(Run3TrksWtofWevTime const& tracks,
                   aod::Collisions const&,
                   aod::BCsWithTimestamps const& bcs)
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

    tofResponse->processSetup(bcs.iteratorAt(0)); // Update the calibration parameters

    for (auto const& pidId : mEnabledParticles) {
      reserveTable(pidId, tracks.size(), false);
    }

    for (auto const& pidId : mEnabledParticlesFull) {
      reserveTable(pidId, tracks.size(), true);
    }

    float resolution = 1.f; // Last resolution assigned
    float nsigma = 0;
    for (auto const& trk : tracks) { // Loop on all tracks
      if (!trk.has_collision()) {    // Track was not assigned, cannot compute NSigma (no event time) -> filling with empty table
        for (auto const& pidId : mEnabledParticles) {
          makeTableEmpty(pidId, false);
        }
        for (auto const& pidId : mEnabledParticlesFull) {
          makeTableEmpty(pidId, true);
        }
        continue;
      }

      for (auto const& pidId : mEnabledParticles) { // Loop on enabled particle hypotheses
        switch (pidId) {
          case kIdxEl: {
            nsigma = responseEl.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDEl);
            break;
          }
          case kIdxMu: {
            nsigma = responseMu.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDMu);
            break;
          }
          case kIdxPi: {
            nsigma = responsePi.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDPi);
            break;
          }
          case kIdxKa: {
            nsigma = responseKa.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDKa);
            break;
          }
          case kIdxPr: {
            nsigma = responsePr.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDPr);
            break;
          }
          case kIdxDe: {
            nsigma = responseDe.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDDe);
            break;
          }
          case kIdxTr: {
            nsigma = responseTr.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDTr);
            break;
          }
          case kIdxHe: {
            nsigma = responseHe.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDHe);
            break;
          }
          case kIdxAl: {
            nsigma = responseAl.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDAl);
            break;
          }
          default:
            LOG(fatal) << "Wrong particle ID for standard tables";
            break;
        }
        if (enableQaHistograms) {
          hnsigma[pidId]->Fill(trk.p(), nsigma);
        }
      }
      for (auto const& pidId : mEnabledParticlesFull) { // Loop on enabled particle hypotheses with full tables
        switch (pidId) {
          case kIdxEl: {
            resolution = responseEl.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseEl.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullEl(resolution, nsigma);
            break;
          }
          case kIdxMu: {
            resolution = responseMu.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseMu.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullMu(resolution, nsigma);
            break;
          }
          case kIdxPi: {
            resolution = responsePi.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responsePi.GetSeparation(tofResponse->parameters, trk);
            tablePIDFullPi(resolution, nsigma);
            break;
          }
          case kIdxKa: {
            resolution = responseKa.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseKa.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullKa(resolution, nsigma);
            break;
          }
          case kIdxPr: {
            resolution = responsePr.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responsePr.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullPr(resolution, nsigma);
            break;
          }
          case kIdxDe: {
            resolution = responseDe.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseDe.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullDe(resolution, nsigma);
            break;
          }
          case kIdxTr: {
            resolution = responseTr.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseTr.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullTr(resolution, nsigma);
            break;
          }
          case kIdxHe: {
            resolution = responseHe.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseHe.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullHe(resolution, nsigma);
            break;
          }
          case kIdxAl: {
            resolution = responseAl.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseAl.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullAl(resolution, nsigma);
            break;
          }
          default:
            LOG(fatal) << "Wrong particle ID for full tables";
            break;
        }
        if (enableQaHistograms) {
          hnsigmaFull[pidId]->Fill(trk.p(), nsigma);
        }
      }
    }
  }
  PROCESS_SWITCH(tofPidMerge, processRun3, "Produce Run 3 Nsigma table. Set to off if the tables are not required, or autoset is on", false);

  template <o2::track::PID::ID pid>
  using ResponseImplementationRun2 = o2::pid::tof::ExpTimes<Run2TrksWtofWevTime::iterator, pid>;
  void processRun2(Run2TrksWtofWevTime const& tracks,
                   aod::Collisions const&,
                   aod::BCsWithTimestamps const& bcs)
  {
    constexpr auto responseEl = ResponseImplementationRun2<PID::Electron>();
    constexpr auto responseMu = ResponseImplementationRun2<PID::Muon>();
    constexpr auto responsePi = ResponseImplementationRun2<PID::Pion>();
    constexpr auto responseKa = ResponseImplementationRun2<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementationRun2<PID::Proton>();
    constexpr auto responseDe = ResponseImplementationRun2<PID::Deuteron>();
    constexpr auto responseTr = ResponseImplementationRun2<PID::Triton>();
    constexpr auto responseHe = ResponseImplementationRun2<PID::Helium3>();
    constexpr auto responseAl = ResponseImplementationRun2<PID::Alpha>();

    tofResponse->processSetup(bcs.iteratorAt(0)); // Update the calibration parameters

    for (auto const& pidId : mEnabledParticles) {
      reserveTable(pidId, tracks.size(), false);
    }

    for (auto const& pidId : mEnabledParticlesFull) {
      reserveTable(pidId, tracks.size(), true);
    }

    float resolution = 1.f; // Last resolution assigned
    float nsigma = 0;
    for (auto const& trk : tracks) { // Loop on all tracks
      if (!trk.has_collision()) {    // Track was not assigned, cannot compute NSigma (no event time) -> filling with empty table
        for (auto const& pidId : mEnabledParticles) {
          makeTableEmpty(pidId, false);
        }
        for (auto const& pidId : mEnabledParticlesFull) {
          makeTableEmpty(pidId, true);
        }
        continue;
      }

      for (auto const& pidId : mEnabledParticles) { // Loop on enabled particle hypotheses
        switch (pidId) {
          case kIdxEl: {
            nsigma = responseEl.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDEl);
            break;
          }
          case kIdxMu: {
            nsigma = responseMu.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDMu);
            break;
          }
          case kIdxPi: {
            nsigma = responsePi.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDPi);
            break;
          }
          case kIdxKa: {
            nsigma = responseKa.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDKa);
            break;
          }
          case kIdxPr: {
            nsigma = responsePr.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDPr);
            break;
          }
          case kIdxDe: {
            nsigma = responseDe.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDDe);
            break;
          }
          case kIdxTr: {
            nsigma = responseTr.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDTr);
            break;
          }
          case kIdxHe: {
            nsigma = responseHe.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDHe);
            break;
          }
          case kIdxAl: {
            nsigma = responseAl.GetSeparation(tofResponse->parameters, trk);
            aod::pidtof_tiny::binning::packInTable(nsigma, tablePIDAl);
            break;
          }
          default:
            LOG(fatal) << "Wrong particle ID for standard tables";
            break;
        }
        if (enableQaHistograms) {
          hnsigma[pidId]->Fill(trk.p(), nsigma);
        }
      }
      for (auto const& pidId : mEnabledParticlesFull) { // Loop on enabled particle hypotheses with full tables
        switch (pidId) {
          case kIdxEl: {
            resolution = responseEl.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseEl.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullEl(resolution, nsigma);
            break;
          }
          case kIdxMu: {
            resolution = responseMu.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseMu.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullMu(resolution, nsigma);
            break;
          }
          case kIdxPi: {
            resolution = responsePi.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responsePi.GetSeparation(tofResponse->parameters, trk);
            tablePIDFullPi(resolution, nsigma);
            break;
          }
          case kIdxKa: {
            resolution = responseKa.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseKa.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullKa(resolution, nsigma);
            break;
          }
          case kIdxPr: {
            resolution = responsePr.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responsePr.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullPr(resolution, nsigma);
            break;
          }
          case kIdxDe: {
            resolution = responseDe.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseDe.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullDe(resolution, nsigma);
            break;
          }
          case kIdxTr: {
            resolution = responseTr.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseTr.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullTr(resolution, nsigma);
            break;
          }
          case kIdxHe: {
            resolution = responseHe.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseHe.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullHe(resolution, nsigma);
            break;
          }
          case kIdxAl: {
            resolution = responseAl.GetExpectedSigma(tofResponse->parameters, trk);
            nsigma = responseAl.GetSeparation(tofResponse->parameters, trk, resolution);
            tablePIDFullAl(resolution, nsigma);
            break;
          }
          default:
            LOG(fatal) << "Wrong particle ID for full tables";
            break;
        }
        if (enableQaHistograms) {
          hnsigmaFull[pidId]->Fill(trk.p(), nsigma);
        }
      }
    }
  }
  PROCESS_SWITCH(tofPidMerge, processRun2, "Produce Run 2 Nsigma table. Set to off if the tables are not required, or autoset is on", false);

  o2::pid::tof::Beta responseBetaRun2;
  void processRun2BetaM(Run2TrksWtofWevTime const& tracks)
  {
    if (!enableTableBeta && !enableTableMass) {
      return;
    }
    float beta = 0.f;
    tablePIDBeta.reserve(tracks.size());
    for (auto const& trk : tracks) {
      beta = responseBetaRun2.GetBeta(trk);
      if (enableTableBeta) {
        tablePIDBeta(beta, responseBetaRun2.GetExpectedSigma(trk));
      }
      if (enableTableMass) {
        if (enableTOFParamsForBetaMass) {
          tablePIDTOFMass(o2::pid::tof::TOFMass::GetTOFMass(trk.tofExpMom() / (1.f + trk.sign() * tofResponse->parameters.getMomentumChargeShift(trk.eta())), beta));
        } else {
          tablePIDTOFMass(o2::pid::tof::TOFMass::GetTOFMass(trk, beta));
        }
      }
    }
  }
  PROCESS_SWITCH(tofPidMerge, processRun2BetaM, "Produce Run 2 Beta and Mass table. Set to off if the tables are not required, or autoset is on", false);

  o2::pid::tof::Beta responseBeta;
  void processRun3BetaM(Run3TrksWtofWevTime const& tracks)
  {
    if (!enableTableBeta && !enableTableMass) {
      return;
    }
    float beta = 0.f;
    tablePIDBeta.reserve(tracks.size());
    for (auto const& trk : tracks) {
      beta = responseBeta.GetBeta(trk);
      if (enableTableBeta) {
        tablePIDBeta(beta,
                     responseBeta.GetExpectedSigma(trk));
      }
      if (enableTableMass) {
        if (enableTOFParamsForBetaMass) {
          tablePIDTOFMass(o2::pid::tof::TOFMass::GetTOFMass(trk.tofExpMom() / (1.f + trk.sign() * tofResponse->parameters.getMomentumChargeShift(trk.eta())), beta));
        } else {
          tablePIDTOFMass(o2::pid::tof::TOFMass::GetTOFMass(trk, beta));
        }
      }
    }
  }
  PROCESS_SWITCH(tofPidMerge, processRun3BetaM, "Produce Run 3 Beta and Mass table. Set to off if the tables are not required, or autoset is on", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(cfgc);
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofSignal>(cfgc)};
  workflow.push_back(adaptAnalysisTask<tofEventTime>(cfgc));
  workflow.push_back(adaptAnalysisTask<tofPidMerge>(cfgc));
  return workflow;
}
