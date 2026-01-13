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
/// \file   mcPidTof.cxx
/// \author Fabrizio Grosa fabrizio.grosa@cern.ch
/// \brief  Task to produce PID tables for TOF split for pi, K, p, de, copied from https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/PID/pidTofMerge.cxx
///         It works only for MC and adds the possibility to apply postcalibrations for MC.
///

#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
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
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <PID/PIDTOF.h>
#include <ReconstructionDataFormats/PID.h>
#include <TOFBase/EventTimeMaker.h>

#include <TGraph.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
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

o2::common::core::MetadataHelper metadataInfo;

// Input data types
using Trks = o2::soa::Join<aod::TracksIU, aod::TracksExtra>;
using Cols = aod::Collisions;
using TrksWtof = soa::Join<Trks, aod::TOFSignal>;
using TrksWtofWevTime = soa::Join<TrksWtof, aod::TOFEvTime, aod::pidEvTimeFlags, aod::McTrackLabels>;

using EvTimeCollisions = soa::Join<Cols, aod::EvSels>;
using EvTimeCollisionsFT0 = soa::Join<EvTimeCollisions, aod::FT0sCorrected>;

// Configuration common to all tasks
struct TOFCalibConfig {
  template <typename CfgType>
  void init(const CfgType& opt)
  {
    mUrl = opt.cfgUrl.value;
    mPathGrpLhcIf = opt.cfgPathGrpLhcIf.value;
    mTimestamp = opt.cfgTimestamp.value;
    mTimeShiftCCDBPathPos = opt.cfgTimeShiftCCDBPathPos.value;
    mTimeShiftCCDBPathNeg = opt.cfgTimeShiftCCDBPathNeg.value;
    mTimeShiftCCDBPathPosMC = opt.cfgTimeShiftCCDBPathPosMC.value;
    mTimeShiftCCDBPathNegMC = opt.cfgTimeShiftCCDBPathNegMC.value;
    mParamFileName = opt.cfgParamFileName.value;
    mParametrizationPath = opt.cfgParametrizationPath.value;
    mReconstructionPass = opt.cfgReconstructionPass.value;
    mReconstructionPassDefault = opt.cfgReconstructionPassDefault.value;
    mFatalOnPassNotAvailable = opt.cfgFatalOnPassNotAvailable.value;
    mEnableTimeDependentResponse = opt.cfgEnableTimeDependentResponse.value;
    mCollisionSystem = opt.cfgCollisionSystem.value;
    mAutoSetProcessFunctions = opt.cfgAutoSetProcessFunctions.value;
  }

  template <typename VType>
  void getCfg(o2::framework::InitContext& initContext, const std::string& name, VType& v, const std::string& task)
  {
    if (!getTaskOptionValue(initContext, task, name, v, false)) {
      LOG(fatal) << "Could not get " << name << " from " << task << " task";
    }
  }

  void inheritFromBaseTask(o2::framework::InitContext& initContext, const std::string& task = "tof-signal")
  {
    mInitMode = 2;
    getCfg(initContext, "ccdb-url", mUrl, task);
    getCfg(initContext, "ccdb-path-grplhcif", mPathGrpLhcIf, task);
    getCfg(initContext, "ccdb-timestamp", mTimestamp, task);
    getCfg(initContext, "timeShiftCCDBPathPos", mTimeShiftCCDBPathPos, task);
    getCfg(initContext, "timeShiftCCDBPathNeg", mTimeShiftCCDBPathNeg, task);
    getCfg(initContext, "timeShiftCCDBPathPosMC", mTimeShiftCCDBPathPosMC, task);
    getCfg(initContext, "timeShiftCCDBPathNegMC", mTimeShiftCCDBPathNegMC, task);
    getCfg(initContext, "paramFileName", mParamFileName, task);
    getCfg(initContext, "parametrizationPath", mParametrizationPath, task);
    getCfg(initContext, "reconstructionPass", mReconstructionPass, task);
    getCfg(initContext, "reconstructionPassDefault", mReconstructionPassDefault, task);
    getCfg(initContext, "fatalOnPassNotAvailable", mFatalOnPassNotAvailable, task);
    getCfg(initContext, "enableTimeDependentResponse", mEnableTimeDependentResponse, task);
    getCfg(initContext, "collisionSystem", mCollisionSystem, task);
    getCfg(initContext, "autoSetProcessFunctions", mAutoSetProcessFunctions, task);
  }
  // @brief Set up the configuration from the calibration object from the init function of the task
  template <typename CCDBObject>
  void initSetup(o2::pid::tof::TOFResoParamsV3& mRespParamsV3,
                 CCDBObject ccdb)
  {
    mInitMode = 1;
    // First we set the CCDB manager
    ccdb->setURL(mUrl);
    ccdb->setTimestamp(mTimestamp);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    // Then the information about the metadata
    if (mReconstructionPass == "metadata") {
      LOG(info) << "Getting pass from metadata";
      if (metadataInfo.isMC()) {
        mReconstructionPass = metadataInfo.get("AnchorPassName");
      } else {
        LOG(fatal) << "This task works only for MC";
      }
      LOG(info) << "Passed autodetect mode for pass. Taking '" << mReconstructionPass << "'";
    }
    LOG(info) << "Using parameter collection, starting from pass '" << mReconstructionPass << "'";

    if (!mParamFileName.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file " << mParamFileName << ", using param: " << mParametrizationPath << " and pass " << mReconstructionPass;
      o2::tof::ParameterCollection paramCollection;
      paramCollection.loadParamFromFile(mParamFileName, mParametrizationPath);
      LOG(info) << "+++ Loaded parameter collection from file +++";
      if (!paramCollection.retrieveParameters(mRespParamsV3, mReconstructionPass)) {
        if (mFatalOnPassNotAvailable) {
          LOG(fatal) << "Pass '" << mReconstructionPass << "' not available in the retrieved object from file";
        } else {
          LOG(warning) << "Pass '" << mReconstructionPass << "' not available in the retrieved object from file, fetching '" << mReconstructionPassDefault << "'";
          if (!paramCollection.retrieveParameters(mRespParamsV3, mReconstructionPassDefault)) {
            paramCollection.print();
            LOG(fatal) << "Cannot get default pass for calibration " << mReconstructionPassDefault;
          } else {
            mRespParamsV3.setResolutionParametrization(paramCollection.getPars(mReconstructionPassDefault));
            mRespParamsV3.setMomentumChargeShiftParameters(paramCollection.getPars(mReconstructionPassDefault));
          }
        }
      } else { // Pass is available, load non standard parameters
        mRespParamsV3.setResolutionParametrization(paramCollection.getPars(mReconstructionPass));
        mRespParamsV3.setMomentumChargeShiftParameters(paramCollection.getPars(mReconstructionPass));
      }
    } else if (!mEnableTimeDependentResponse) { // Loading it from CCDB
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << mParametrizationPath << " for timestamp " << mTimestamp;
      o2::tof::ParameterCollection* paramCollection = ccdb->template getForTimeStamp<o2::tof::ParameterCollection>(mParametrizationPath, mTimestamp);
      if (!paramCollection->retrieveParameters(mRespParamsV3, mReconstructionPass)) { // Attempt at loading the parameters with the pass defined
        if (mFatalOnPassNotAvailable) {
          LOG(fatal) << "Pass '" << mReconstructionPass << "' not available in the retrieved CCDB object";
        } else {
          LOG(warning) << "Pass '" << mReconstructionPass << "' not available in the retrieved CCDB object, fetching '" << mReconstructionPassDefault << "'";
          if (!paramCollection->retrieveParameters(mRespParamsV3, mReconstructionPassDefault)) {
            paramCollection->print();
            LOG(fatal) << "Cannot get default pass for calibration " << mReconstructionPassDefault;
          } else {
            mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPassDefault));
            mRespParamsV3.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPassDefault));
          }
        }
      } else { // Pass is available, load non standard parameters
        mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPass));
        mRespParamsV3.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPass));
      }
    }

    // Loading additional calibration objects
    std::map<std::string, std::string> metadata;
    if (!mReconstructionPass.empty()) {
      metadata["RecoPassName"] = mReconstructionPass;
    }
    auto updateTimeShift = [&](const std::string& nameShift, bool isPositive) {
      if (nameShift.empty()) {
        return;
      }
      const bool isFromFile = nameShift.find(".root") != std::string::npos;
      if (isFromFile) {
        LOG(info) << "Initializing the time shift for " << (isPositive ? "positive" : "negative") << " from file '" << nameShift << "'";
        mRespParamsV3.setTimeShiftParameters(nameShift, "ccdb_object", isPositive);
      } else if (!mEnableTimeDependentResponse) { // If the response is fixed fetch it at the init time
        LOG(info) << "Initializing the time shift for " << (isPositive ? "positive" : "negative")
                  << " from ccdb '" << nameShift << "' and timestamp " << mTimestamp
                  << " and pass '" << mReconstructionPass << "'";
        mRespParamsV3.setTimeShiftParameters(ccdb->template getSpecific<TGraph>(nameShift, mTimestamp, metadata), isPositive);
      }
      LOG(info) << " test getTimeShift at 0 " << (isPositive ? "pos" : "neg") << ": "
                << mRespParamsV3.getTimeShift(0, isPositive);
    };

    const std::string nameShiftPos = metadataInfo.isMC() ? mTimeShiftCCDBPathPosMC : mTimeShiftCCDBPathPos;
    updateTimeShift(nameShiftPos, true);
    const std::string nameShiftNeg = metadataInfo.isMC() ? mTimeShiftCCDBPathNegMC : mTimeShiftCCDBPathNeg;
    updateTimeShift(nameShiftNeg, true);

    // Calibration object is defined
    LOG(info) << "Parametrization at init time:";
    mRespParamsV3.printFullConfig();
  }

  template <typename CCDBObject, typename BcType>
  void processSetup(o2::pid::tof::TOFResoParamsV3& mRespParamsV3,
                    CCDBObject ccdb,
                    const BcType& bc)
  {
    LOG(debug) << "Processing setup for run number " << bc.runNumber() << " from run " << mLastRunNumber;
    // First we check if this run number was already processed
    if (mLastRunNumber == bc.runNumber()) {
      return;
    }
    mLastRunNumber = bc.runNumber();
    mTimestamp = bc.timestamp();

    // Check the beam type
    if (mCollisionSystem == -1) {
      o2::parameters::GRPLHCIFData* grpo = ccdb->template getSpecific<o2::parameters::GRPLHCIFData>(mPathGrpLhcIf,
                                                                                                    mTimestamp);
      mCollisionSystem = CollisionSystemType::getCollisionTypeFromGrp(grpo);
    } else {
      LOG(debug) << "Not setting collisions system as already set to " << mCollisionSystem << " " << CollisionSystemType::getCollisionSystemName(mCollisionSystem);
    }

    if (!mEnableTimeDependentResponse) {
      return;
    }
    LOG(info) << "Updating parametrization from path '" << mParametrizationPath << "' and timestamp " << mTimestamp << " and reconstruction pass '" << mReconstructionPass << "' for run number " << bc.runNumber();
    if (mParamFileName.empty()) { // Not loading if parametrization was taken from file
      LOG(info) << "Updating parametrization from ccdb";
      const o2::tof::ParameterCollection* paramCollection = ccdb->template getSpecific<o2::tof::ParameterCollection>(mParametrizationPath, mTimestamp);
      if (!paramCollection->retrieveParameters(mRespParamsV3, mReconstructionPass)) {
        if (mFatalOnPassNotAvailable) {
          LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", mReconstructionPass.data());
        } else {
          LOGF(warning, "Pass '%s' not available in the retrieved CCDB object, fetching '%s'", mReconstructionPass.data(), mReconstructionPassDefault.data());
          if (!paramCollection->retrieveParameters(mRespParamsV3, mReconstructionPassDefault)) {
            paramCollection->print();
            LOG(fatal) << "Cannot get default pass for calibration " << mReconstructionPassDefault;
          } else { // Found the default case
            mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPassDefault));
            mRespParamsV3.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPassDefault));
          }
        }
      } else { // Found the non default case
        mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPass));
        mRespParamsV3.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPass));
      }
    }

    // Loading additional calibration objects
    std::map<std::string, std::string> metadata;
    if (!mReconstructionPass.empty()) {
      metadata["RecoPassName"] = mReconstructionPass;
    }
    auto updateTimeShift = [&](const std::string& nameShift, bool isPositive) {
      if (nameShift.empty()) {
        return;
      }
      const bool isFromFile = nameShift.find(".root") != std::string::npos;
      if (isFromFile) {
        return;
      }
      LOG(info) << "Updating the time shift for " << (isPositive ? "positive" : "negative")
                << " from ccdb '" << nameShift << "' and timestamp " << mTimestamp
                << " and pass '" << mReconstructionPass << "'";
      mRespParamsV3.setTimeShiftParameters(ccdb->template getSpecific<TGraph>(nameShift, mTimestamp, metadata), isPositive);
      LOG(info) << " test getTimeShift at 0 " << (isPositive ? "pos" : "neg") << ": "
                << mRespParamsV3.getTimeShift(0, isPositive);
    };

    updateTimeShift(metadataInfo.isMC() ? mTimeShiftCCDBPathPosMC : mTimeShiftCCDBPathPos, true);
    updateTimeShift(metadataInfo.isMC() ? mTimeShiftCCDBPathNegMC : mTimeShiftCCDBPathNeg, false);

    LOG(info) << "Parametrization at setup time:";
    mRespParamsV3.printFullConfig();
  }

  bool autoSetProcessFunctions() const { return mAutoSetProcessFunctions; }
  int collisionSystem() const { return mCollisionSystem; }

 private:
  int mLastRunNumber = -1; // Last run number for which the calibration was loaded
  int mInitMode = 0;       // 0: no init, 1: init, 2: inherit

  // Configurable options
  std::string mUrl;
  std::string mPathGrpLhcIf;
  int64_t mTimestamp{0};
  std::string mTimeShiftCCDBPathPos;
  std::string mTimeShiftCCDBPathNeg;
  std::string mTimeShiftCCDBPathPosMC;
  std::string mTimeShiftCCDBPathNegMC;
  std::string mParamFileName;
  std::string mParametrizationPath;
  std::string mReconstructionPass;
  std::string mReconstructionPassDefault;
  bool mFatalOnPassNotAvailable{false};
  bool mEnableTimeDependentResponse{false};
  int mCollisionSystem{-1};
  bool mAutoSetProcessFunctions{false};
};

// Part 1 TOF signal definition

/// Selection criteria for tracks used for TOF event time
bool isTrackGoodMatchForTOFPID(const Trks::iterator& tr)
{
  return tr.hasTOF();
}

/// Task to produce the TOF signal from the trackTime information
struct TofSignal {
  // Tables to produce
  o2::framework::Produces<o2::aod::TOFSignal> table;
  o2::framework::Produces<o2::aod::pidTOFFlags> tableFlags;
  // Running flags
  bool enableTableTOFSignal = false;   // Flag to check if the TOF signal table is requested or not
  bool enableTablepidTOFFlags = false; // Flag to check if the TOF signal flags table is requested or not
  // Output histograms
  Configurable<bool> enableQaHistograms{"enableQaHistograms", false, "Flag to enable the QA histograms"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // Detector response and input parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
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

  TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration

  void init(o2::framework::InitContext& initContext)
  {
    mTOFCalibConfig.init(cfg);
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
    if (!enableTableTOFSignal && !enableTablepidTOFFlags) {
      LOG(info) << "No table or process is enabled. Disabling task";
      return;
    }

    mTOFCalibConfig.initSetup(mRespParamsV3, ccdb); // Getting the parametrization parameters
    if (!enableQaHistograms) {
      return;
    }
    histos.add("tofSignal", "tofSignal", kTH1D, {{1000, -1000, 1000000, "tofSignal (ps)"}});
    if (enableTablepidTOFFlags) {
      histos.add("goodForPIDFlags", "goodForPIDFlags", kTH1D, {{3, 0, 3, "flags"}});
    }
  }

  void process(Trks const& tracks)
  {
    if (!enableTableTOFSignal) {
      return;
    }
    table.reserve(tracks.size());
    if (enableTablepidTOFFlags) {
      tableFlags.reserve(tracks.size());
    }
    for (auto const& trk : tracks) {
      const float& sig = o2::pid::tof::TOFSignal<Trks::iterator>::GetTOFSignal(trk);
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
};

/// Selection criteria for tracks used for TOF event time
float trackSampleMinMomentum = 0.5f;
float trackSampleMaxMomentum = 2.f;
template <typename TrackType>
bool filterForTOFEventTime(const TrackType& tr)
{
  return (tr.hasTOF() &&
          tr.p() > trackSampleMinMomentum && tr.p() < trackSampleMaxMomentum &&
          tr.hasITS() &&
          tr.hasTPC() &&
          (tr.trackType() == o2::aod::track::TrackTypeEnum::Track || tr.trackType() == o2::aod::track::TrackTypeEnum::TrackIU));
} // accept all

/// Specialization of TOF event time maker
template <typename TrackType,
          bool (*TrackFilter)(const TrackType&),
          template <typename T, o2::track::PID::ID> typename Response,
          typename TrackTypeContainer,
          typename ResponseParametersType>
o2::tof::eventTimeContainer evTimeMakerForTracks(const TrackTypeContainer& tracks,
                                                 const ResponseParametersType& responseParameters,
                                                 const float& diamond = 6.0)
{
  return o2::tof::evTimeMakerFromParam<TrackTypeContainer, TrackType, TrackFilter, Response, ResponseParametersType>(tracks, responseParameters, diamond);
}

// Part 2 event time definition

/// Task to produce the TOF event time table
struct TofEventTime {
  // Tables to produce
  Produces<o2::aod::TOFEvTime> tableEvTime;
  Produces<o2::aod::EvTimeTOFOnly> tableEvTimeTOFOnly;
  Produces<o2::aod::pidEvTimeFlags> tableFlags;
  static constexpr bool RemoveTofEvTimeBias = true; // Flag to subtract the Ev. Time bias for low multiplicity events with TOF
  static constexpr float Diamond = 6.0;             // Collision diamond used in the estimation of the TOF event time
  static constexpr float ErrDiamond = Diamond * 33.356409f;
  static constexpr float WeightDiamond = 1.f / (ErrDiamond * ErrDiamond);

  bool enableTableTOFEvTime = false;
  bool enableTableEvTimeTOFOnly = false;
  // Detector response and input parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration

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
    mTOFCalibConfig.inheritFromBaseTask(initContext);
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

    if (metadataInfo.isFullyDefined()) {
      if (!metadataInfo.isRun3()) {
        LOG(fatal) << "Metadata says it is Run2, but this task supports only Run3";
      }
    }

    trackSampleMinMomentum = minMomentum;
    trackSampleMaxMomentum = maxMomentum;
    LOG(info) << "Configuring track sample for TOF ev. time: " << trackSampleMinMomentum << " < p < " << trackSampleMaxMomentum;

    if (sel8TOFEvTime.value) {
      LOG(info) << "TOF event time will be computed for collisions that pass the event selection only!";
    }
    mTOFCalibConfig.initSetup(mRespParamsV3, ccdb); // Getting the parametrization parameters

    o2::tof::eventTimeContainer::setMaxNtracksInSet(maxNtracksInSet.value);
    o2::tof::eventTimeContainer::printConfig();
  }

  ///
  /// Process function to prepare the event for each track on Run 3 data without the FT0
  // Define slice per collision
  Preslice<TrksWtof> perCollision = aod::track::collisionId;
  template <o2::track::PID::ID Pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksWtof::iterator, Pid>;
  void process(TrksWtof const& tracks,
               aod::FT0s const&,
               EvTimeCollisionsFT0 const&,
               aod::BCsWithTimestamps const& bcs)
  {
    if (!enableTableTOFEvTime) {
      return;
    }
    LOG(debug) << "Processing data for TOF event time";

    tableEvTime.reserve(tracks.size());
    tableFlags.reserve(tracks.size());
    if (enableTableEvTimeTOFOnly) {
      tableEvTimeTOFOnly.reserve(tracks.size());
    }

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
      int lastCollisionId = -1;                                                                               // Last collision ID analysed
      for (auto const& t : tracks) {                                                                          // Loop on collisions
        if (!t.has_collision() || ((sel8TOFEvTime.value) && !t.collision_as<EvTimeCollisionsFT0>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
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
        const auto evTimeMakerTOF = evTimeMakerForTracks<TrksWtof::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV3, Diamond);

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
          if constexpr (RemoveTofEvTimeBias) {
            evTimeMakerTOF.removeBias<TrksWtof::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, t0TOF[0], t0TOF[1], 2);
          }
          if (t0TOF[1] < ErrDiamond && (maxEvTimeTOF <= 0 || std::abs(t0TOF[0]) < maxEvTimeTOF)) {
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

          if (sumOfWeights < WeightDiamond) { // avoiding sumOfWeights = 0 or worse that diamond
            eventTime = 0;
            sumOfWeights = WeightDiamond;
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
      int lastCollisionId = -1;                                                                            // Last collision ID analysed
      for (auto const& t : tracks) {                                                                       // Loop on collisions
        if (!t.has_collision() || ((sel8TOFEvTime.value) && !t.collision_as<EvTimeCollisions>().sel8())) { // Track was not assigned, cannot compute event time or event did not pass the event selection
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
        const auto evTimeMakerTOF = evTimeMakerForTracks<TrksWtof::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, mRespParamsV3, Diamond);
        int nGoodTracksForTOF = 0;
        float et = evTimeMakerTOF.mEventTime;
        float erret = evTimeMakerTOF.mEventTimeError;

        for (auto const& trk : tracksInCollision) { // Loop on Tracks
          if constexpr (RemoveTofEvTimeBias) {
            evTimeMakerTOF.removeBias<TrksWtof::iterator, filterForTOFEventTime>(trk, nGoodTracksForTOF, et, erret, 2);
          }
          uint8_t flags = 0;
          if (erret < ErrDiamond && (maxEvTimeTOF <= 0.f || std::abs(et) < maxEvTimeTOF)) {
            flags |= o2::aod::pidflags::enums::PIDFlags::EvTimeTOF;
          } else {
            et = 0.f;
            erret = ErrDiamond;
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
};

// Part 3 Nsigma computation

static constexpr int IdxPi = 2;
static constexpr int IdxKa = 3;
static constexpr int IdxPr = 4;
static constexpr int IdxDe = 5;

/// Task to produce the response table
struct McPidTof {
  // Tables to produce
  Produces<o2::aod::pidTOFPi> tablePIDPi;
  Produces<o2::aod::pidTOFKa> tablePIDKa;
  Produces<o2::aod::pidTOFPr> tablePIDPr;
  Produces<o2::aod::pidTOFDe> tablePIDDe;

  // Tables to produce (full)
  Produces<o2::aod::pidTOFFullPi> tablePIDFullPi;
  Produces<o2::aod::pidTOFFullKa> tablePIDFullKa;
  Produces<o2::aod::pidTOFFullPr> tablePIDFullPr;
  Produces<o2::aod::pidTOFFullDe> tablePIDFullDe;

  // Detector response parameters
  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  TOFCalibConfig mTOFCalibConfig; // TOF Calib configuration
  Configurable<bool> enableQaHistograms{"enableQaHistograms", false, "Flag to enable the QA histograms"};

  // Histograms for QA
  std::array<std::shared_ptr<TH2>, nSpecies> hnSigma;
  std::array<std::shared_ptr<TH2>, nSpecies> hnSigmaFull;

  // postcalibrations to overcome MC FT0 timing issue
  std::map<int, TGraph*> gMcPostCalibMean;
  std::map<int, TGraph*> gMcPostCalibSigma;
  int currentRun{0};
  struct : ConfigurableGroup {
    std::string prefix = "mcRecalib";
    Configurable<bool> enable{"enable", false, "enable MC recalibration for Pi/Ka/Pr"};
    Configurable<std::string> ccdbPath{"ccdbPath", "Users/f/fgrosa/RecalibmcPidTof/", "path for MC recalibration objects in CCDB"};
  } mcRecalib;

  // list of productions for which the postcalibrations is needed (bug in FT0 digitisation)
  const std::vector<std::string> prodPostCalib = {"LHC24d3a", "LHC24d3b", "LHC24e3", "LHC24g5", "LHC24g6", "LHC24h2", "LHC24i1", "LHC24i2", "LHC24i3", "LHC24i4", "LHC24j6", "LHC24k3", "LHC24f3c", "LHC25e2", "LHC24h1b", "LHC25e4", "LHC25f5", "LHC25e8", "LHC25e9", "LHC25e10", "LHC25e11", "LHC23k4"};
  bool enableMcRecalib{false};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Running variables
  std::vector<int> mEnabledParticles;     // Vector of enabled PID hypotheses to loop on when making tables
  std::vector<int> mEnabledParticlesFull; // Vector of enabled PID hypotheses to loop on when making full tables
  void init(o2::framework::InitContext& initContext)
  {
    mTOFCalibConfig.inheritFromBaseTask(initContext);
    // Checking the tables are requested in the workflow and enabling them (only pi, K, p)
    std::array<int, 4> supportedSpecies = {IdxPi, IdxKa, IdxPr, IdxDe};
    for (auto iSpecie{0u}; iSpecie < supportedSpecies.size(); ++iSpecie) {
      // First checking tiny
      int flag = -1;
      enableFlagIfTableRequired(initContext, "pidTOF" + particleNames[supportedSpecies[iSpecie]], flag);
      if (flag == 1) {
        mEnabledParticles.push_back(supportedSpecies[iSpecie]);
      }

      // Then check full
      flag = -1;
      enableFlagIfTableRequired(initContext, "pidTOFFull" + particleNames[supportedSpecies[iSpecie]], flag);
      if (flag == 1) {
        mEnabledParticlesFull.push_back(supportedSpecies[iSpecie]);
      }
    }
    if (mEnabledParticlesFull.empty() && mEnabledParticles.empty()) {
      LOG(info) << "No PID tables are required, disabling process function";
      doprocessFillTables.value = false;
      doprocessDummy.value = true;
      return;
    }
    if (metadataInfo.isFullyDefined()) {
      if (!metadataInfo.isRun3()) {
        LOG(fatal) << "Metadata says it is Run2, but this task supports only Run3 data";
      }
    }
    mTOFCalibConfig.initSetup(mRespParamsV3, ccdb); // Getting the parametrization parameters

    // Printing enabled tables and enabling QA histograms if needed
    LOG(info) << "++ Enabled tables:";
    const AxisSpec pAxis{100, 0, 5, "#it{p} (GeV/#it{c})"};
    const AxisSpec nSigmaAxis{100, -10, 10, "N_{#sigma}^{TOF}"};
    for (const int& iSpecie : mEnabledParticles) {
      LOG(info) << "++  pidTOF" << particleNames[iSpecie] << " is enabled";
      if (!enableQaHistograms) {
        continue;
      }
      hnSigma[iSpecie] = histos.add<TH2>(Form("nSigma/%s", particleNames[iSpecie].c_str()), Form("N_{#sigma}^{TOF}(%s)", particleNames[iSpecie].c_str()), kTH2F, {pAxis, nSigmaAxis});
    }
    for (const int& iSpecie : mEnabledParticlesFull) {
      LOG(info) << "++  pidTOFFull" << particleNames[iSpecie] << " is enabled";
      if (!enableQaHistograms) {
        continue;
      }
      hnSigmaFull[iSpecie] = histos.add<TH2>(Form("nSigmaFull/%s", particleNames[iSpecie].c_str()), Form("N_{#sigma}^{TOF}(%s)", particleNames[iSpecie].c_str()), kTH2F, {pAxis, nSigmaAxis});
    }

    enableMcRecalib = mcRecalib.enable;
  }

  // Reserves an empty table for the given particle ID with size of the given track table
  void reserveTable(const int id, const int64_t& size, const bool fullTable = false)
  {
    switch (id) {
      case IdxPi: {
        if (fullTable) {
          tablePIDFullPi.reserve(size);
        } else {
          tablePIDPi.reserve(size);
        }
        break;
      }
      case IdxKa: {
        if (fullTable) {
          tablePIDFullKa.reserve(size);
        } else {
          tablePIDKa.reserve(size);
        }
        break;
      }
      case IdxPr: {
        if (fullTable) {
          tablePIDFullPr.reserve(size);
        } else {
          tablePIDPr.reserve(size);
        }
        break;
      }
      case IdxDe: {
        if (fullTable) {
          tablePIDFullDe.reserve(size);
        } else {
          tablePIDDe.reserve(size);
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
      case IdxPi:
        if (fullTable) {
          tablePIDFullPi(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDPi);
        }
        break;
      case IdxKa:
        if (fullTable) {
          tablePIDFullKa(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDKa);
        }
        break;
      case IdxPr:
        if (fullTable) {
          tablePIDFullPr(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDPr);
        }
        break;
      case IdxDe:
        if (fullTable) {
          tablePIDFullDe(-999.f, -999.f);
        } else {
          aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDDe);
        }
        break;
      default:
        LOG(fatal) << "Wrong particle ID in makeTableEmpty() for " << (fullTable ? "full" : "tiny") << " tables";
        break;
    }
  }

  /// Retrieve MC postcalibration objects from CCDB
  /// \param timestamp timestamp
  void retrieveMcPostCalibFromCcdb(int64_t timestamp)
  {
    std::map<std::string, std::string> metadata;
    if (metadataInfo.isFullyDefined()) {
      metadata["RecoPassName"] = metadataInfo.get("AnchorPassName");
      if (std::find(prodPostCalib.begin(), prodPostCalib.end(), metadataInfo.get("LPMProductionTag")) == prodPostCalib.end()) {
        enableMcRecalib = false;
        LOGP(warn, "Nsigma postcalibrations turned off for {} (new MC productions have FT0 digitisation fixed)", metadataInfo.get("LPMProductionTag"));
        return;
      }
    } else {
      LOGP(error, "Impossible to read metadata! Using default calibrations (2022 apass7)");
      metadata["RecoPassName"] = "";
    }
    auto* calibList = ccdb->getSpecific<TList>(mcRecalib.ccdbPath, timestamp, metadata);
    std::vector<int> updatedSpecies{};
    for (auto const& pidId : mEnabledParticles) { // Loop on enabled particle hypotheses (tiny)
      gMcPostCalibMean[pidId] = reinterpret_cast<TGraph*>(calibList->FindObject(Form("Mean%s", particleNames[pidId].data())));
      gMcPostCalibSigma[pidId] = reinterpret_cast<TGraph*>(calibList->FindObject(Form("Sigma%s", particleNames[pidId].data())));
      updatedSpecies.push_back(pidId);
    }
    for (auto const& pidId : mEnabledParticlesFull) { // Loop on enabled particle hypotheses (full)
      if (std::find(updatedSpecies.begin(), updatedSpecies.end(), pidId) != updatedSpecies.end()) {
        continue;
      }
      gMcPostCalibMean[pidId] = reinterpret_cast<TGraph*>(calibList->FindObject(Form("Mean%s", particleNames[pidId].data())));
      gMcPostCalibSigma[pidId] = reinterpret_cast<TGraph*>(calibList->FindObject(Form("Sigma%s", particleNames[pidId].data())));
    }
  }

  /// Apply MC postcalibrations
  /// \param pidId particle id
  /// \param pt track pT
  template <typename T>
  T applyMcRecalib(int pidId, T trackPt, T nSigma)
  {
    if (nSigma < -998) {
      return nSigma;
    }

    float shift{0.f}, scaleWidth{0.f};
    int const nPoints = gMcPostCalibMean[pidId]->GetN();
    double const ptMin = gMcPostCalibMean[pidId]->GetX()[0];
    double const ptMax = gMcPostCalibMean[pidId]->GetX()[nPoints - 1];
    if (trackPt < ptMin) {
      shift = gMcPostCalibMean[pidId]->Eval(ptMin);
      scaleWidth = gMcPostCalibSigma[pidId]->Eval(ptMin);
    } else if (trackPt > ptMax) {
      shift = gMcPostCalibMean[pidId]->Eval(ptMax);
      scaleWidth = gMcPostCalibSigma[pidId]->Eval(ptMax);
    } else {
      shift = gMcPostCalibMean[pidId]->Eval(trackPt);
      scaleWidth = gMcPostCalibSigma[pidId]->Eval(trackPt);
    }

    T nSigmaCorr = (nSigma - shift) / scaleWidth;
    return nSigmaCorr;
  }

  void processDummy(Trks const&) {}
  PROCESS_SWITCH(McPidTof, processDummy, "Dummy process function", false);

  template <o2::track::PID::ID Pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<TrksWtofWevTime::iterator, Pid>;
  void processFillTables(TrksWtofWevTime const& tracks,
                         Cols const&,
                         aod::BCsWithTimestamps const& bcs,
                         aod::McParticles const&)
  {
    constexpr auto ResponsePi = ResponseImplementation<PID::Pion>();
    constexpr auto ResponseKa = ResponseImplementation<PID::Kaon>();
    constexpr auto ResponsePr = ResponseImplementation<PID::Proton>();
    constexpr auto ResponseDe = ResponseImplementation<PID::Deuteron>();

    mTOFCalibConfig.processSetup(mRespParamsV3, ccdb, bcs.iteratorAt(0)); // Update the calibration parameters

    for (auto const& pidId : mEnabledParticles) {
      reserveTable(pidId, tracks.size(), false);
    }

    for (auto const& pidId : mEnabledParticlesFull) {
      reserveTable(pidId, tracks.size(), true);
    }

    float resolution = 1.f; // Last resolution assigned
    float nSigma = 0;
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

      if (enableMcRecalib) {
        auto runNumber = trk.collision().bc_as<aod::BCsWithTimestamps>().runNumber();
        if (runNumber != currentRun) {
          // update postcalibration files
          auto timestamp = trk.collision().bc_as<aod::BCsWithTimestamps>().timestamp();
          retrieveMcPostCalibFromCcdb(timestamp);
        }
        currentRun = runNumber;
      }

      for (auto const& pidId : mEnabledParticles) { // Loop on enabled particle hypotheses
        switch (pidId) {
          case IdxPi: {
            nSigma = ResponsePi.GetSeparation(mRespParamsV3, trk);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == kPiPlus) { // we rescale only true signal
                nSigma = applyMcRecalib(pidId, trk.pt(), nSigma);
              }
            }
            aod::pidtof_tiny::binning::packInTable(nSigma, tablePIDPi);
            break;
          }
          case IdxKa: {
            nSigma = ResponseKa.GetSeparation(mRespParamsV3, trk);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == kKPlus) { // we rescale only true signal
                nSigma = applyMcRecalib(pidId, trk.pt(), nSigma);
              }
            }
            aod::pidtof_tiny::binning::packInTable(nSigma, tablePIDKa);
            break;
          }
          case IdxPr: {
            nSigma = ResponsePr.GetSeparation(mRespParamsV3, trk);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == kProton) { // we rescale only true signal
                nSigma = applyMcRecalib(pidId, trk.pt(), nSigma);
              }
            }
            aod::pidtof_tiny::binning::packInTable(nSigma, tablePIDPr);
            break;
          }
          case IdxDe: {
            nSigma = ResponseDe.GetSeparation(mRespParamsV3, trk);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == o2::constants::physics::kDeuteron) { // we rescale only true signal
                nSigma = applyMcRecalib(IdxPr, trk.pt(), nSigma);                              // FIXME: currently postcalibrations for protons applied to deuterons, to be checked
              }
            }
            aod::pidtof_tiny::binning::packInTable(nSigma, tablePIDDe);
            break;
          }
          default:
            LOG(fatal) << "Wrong particle ID for standard tables";
            break;
        }
        if (enableQaHistograms) {
          hnSigma[pidId]->Fill(trk.p(), nSigma);
        }
      }

      for (auto const& pidId : mEnabledParticlesFull) { // Loop on enabled particle hypotheses with full tables
        switch (pidId) {
          case IdxPi: {
            resolution = ResponsePi.GetExpectedSigma(mRespParamsV3, trk);
            nSigma = ResponsePi.GetSeparation(mRespParamsV3, trk);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == kPiPlus) { // we rescale only true signal
                nSigma = applyMcRecalib(pidId, trk.pt(), nSigma);
              }
            }
            tablePIDFullPi(resolution, nSigma);
            break;
          }
          case IdxKa: {
            resolution = ResponseKa.GetExpectedSigma(mRespParamsV3, trk);
            nSigma = ResponseKa.GetSeparation(mRespParamsV3, trk, resolution);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == kKPlus) { // we rescale only true signal
                nSigma = applyMcRecalib(pidId, trk.pt(), nSigma);
              }
            }
            tablePIDFullKa(resolution, nSigma);
            break;
          }
          case IdxPr: {
            resolution = ResponsePr.GetExpectedSigma(mRespParamsV3, trk);
            nSigma = ResponsePr.GetSeparation(mRespParamsV3, trk, resolution);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == kProton) { // we rescale only true signal
                nSigma = applyMcRecalib(pidId, trk.pt(), nSigma);
              }
            }
            tablePIDFullPr(resolution, nSigma);
            break;
          }
          case IdxDe: {
            resolution = ResponseDe.GetExpectedSigma(mRespParamsV3, trk);
            nSigma = ResponseDe.GetSeparation(mRespParamsV3, trk, resolution);
            if (enableMcRecalib && trk.has_mcParticle()) {
              if (std::abs(trk.mcParticle().pdgCode()) == o2::constants::physics::kDeuteron) { // we rescale only true signal
                nSigma = applyMcRecalib(IdxPr, trk.pt(), nSigma);                              // FIXME: currently postcalibrations for protons applied to deuterons, to be checked
              }
            }
            tablePIDFullDe(resolution, nSigma);
            break;
          }
          default:
            LOG(fatal) << "Wrong particle ID for full tables";
            break;
        }
        if (enableQaHistograms) {
          hnSigmaFull[pidId]->Fill(trk.p(), nSigma);
        }
      }
    }
  }
  PROCESS_SWITCH(McPidTof, processFillTables, "Process with table filling", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);
  auto workflow = WorkflowSpec{adaptAnalysisTask<TofSignal>(cfgc)};
  workflow.push_back(adaptAnalysisTask<TofEventTime>(cfgc));
  workflow.push_back(adaptAnalysisTask<McPidTof>(cfgc));
  return workflow;
}
