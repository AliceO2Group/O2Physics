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
/// \file   PIDTOFParamService.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  30/06/2025
/// \brief  Implementation of the TOF PID service for the detector response
///

#ifndef COMMON_CORE_PID_PIDTOFPARAMSERVICE_H_
#define COMMON_CORE_PID_PIDTOFPARAMSERVICE_H_

#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/TableHelper.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsTOF/ParameterContainers.h"
#include "Framework/DataTypes.h"
#include "Framework/PID.h"
#include "Framework/Plugins.h"
#include "ReconstructionDataFormats/PID.h"

#include <string>

namespace o2::pid::tof
{

struct TOFResponseImpl {
  static o2::pid::tof::TOFResoParamsV3 parameters;      // TOF response parameters for the expected resolution
  static o2::common::core::MetadataHelper metadataInfo; // Metadata information used for the TOF response

  /// Initialize the TOF response parameters in the init function of each task
  /// \param ccdb Pointer to the CCDB manager
  /// \param initContext Initialization context
  /// \note This function should be called in the init function of each task that uses the TOF response
  /// \note The parameters are loaded from the CCDB and stored in the static variable `parameters`
  /// \note The metadata information is also initialized in this function
  void initSetup(o2::ccdb::BasicCCDBManager* ccdb, o2::framework::InitContext& initContext);

  /// Initialize the TOF response parameters in the init function of each task
  /// \param ccdb Service pointer to the CCDB manager
  template <typename T>
  void initSetup(T ccdb, o2::framework::InitContext& initContext)
  {
    initSetup(ccdb.operator->(), initContext);
  }

  /// Initialize the TOF response parameters in the process function of each task, should be called only at least once per run
  /// \param runNumber Run number for which the calibration is loaded
  /// \param timeStamp Timestamp for which the calibration is loaded
  /// \note This function should be called in the process function of each task that uses the TOF response
  /// \note The parameters are loaded from the CCDB and stored in the static variable `parameters`
  /// \note The metadata information is also initialized in this function
  void processSetup(const int runNumber, const int64_t timeStamp);

  /// Initialize the TOF response parameters in the process function of each task, should be called only at least once per run
  /// \param bc Bunch crossing containing the run number and timestamp for which the calibration is loaded
  template <typename T>
  void processSetup(const T& bc)
  {
    processSetup(bc.runNumber(), bc.timestamp());
  }

  template <o2::track::PID::ID id>
  static float expectedSigma(const float tofSignal,
                             const float tofExpMom,
                             const float momentum,
                             const float eta,
                             const float tofEvTimeErr,
                             const o2::pid::tof::TOFResoParamsV3& params = parameters)
  {
    if (!mIsInit) {
      LOG(fatal) << "TOF response parameters not initialized, call initSetup() first";
    }
    if (mLastRunNumber < 0) {
      LOG(fatal) << "TOF response parameters not initialized, call processSetup() first";
    }
    if (tofSignal <= 0.f) {
      // return o2::pid::tof::defaultReturnValue;
    }
    if (tofExpMom <= 0.f) {
      return o2::pid::tof::defaultReturnValue;
    }
    if (momentum <= 0) {
      return o2::pid::tof::defaultReturnValue;
    }
    const float trackingReso = params.getResolution<id>(momentum, eta);
    const float tofReso = params.getParameter(4);
    if (trackingReso > 0) {
      return std::sqrt(trackingReso * trackingReso +
                       tofReso * tofReso +
                       tofEvTimeErr * tofEvTimeErr);
    }
    constexpr float MassSquared = o2::track::pid_constants::sMasses2[id];
    const float dpp = params.getParameter(0) +
                      params.getParameter(1) * momentum +
                      params.getParameter(2) * o2::constants::physics::MassElectron / momentum;
    const float sigma = dpp * tofSignal / (1. + momentum * momentum / (MassSquared));
    return std::sqrt(sigma * sigma +
                     params.getParameter(3) * params.getParameter(3) / momentum / momentum +
                     tofReso * tofReso +
                     tofEvTimeErr * tofEvTimeErr);
  }

  template <o2::track::PID::ID id, typename TrackType>
  static float expectedSigma(const TrackType& track, const o2::pid::tof::TOFResoParamsV3& params = parameters)
  {
    return expectedSigma<id>(track.tofSignal(), track.tofExpMom(), track.p(), track.eta(), track.tofEvTimeErr(), params);
  }

  template <o2::track::PID::ID id>
  static float nSigma(const float tofSignal,
                      const float tofExpMom,
                      const float length,
                      const float momentum,
                      const float eta,
                      const float tofEvTime,
                      const float tofEvTimeErr,
                      const o2::pid::tof::TOFResoParamsV3& params = parameters)
  {
    if (tofSignal <= 0.f) {
      return o2::pid::tof::defaultReturnValue;
    }
    if (tofExpMom <= 0.f) {
      return o2::pid::tof::defaultReturnValue;
    }
    if (momentum <= 0) {
      return o2::pid::tof::defaultReturnValue;
    }

    const float resolution = expectedSigma<id>(tofSignal, tofExpMom, momentum, eta, tofEvTimeErr, params);
    const float expTime = o2::framework::pid::tof::MassToExpTime(tofExpMom,
                                                                 length,
                                                                 o2::track::pid_constants::sMasses2[id]);
    const float delta = tofSignal - tofEvTime - expTime;
    return delta / resolution;
  }

  template <o2::track::PID::ID id, typename TrackType>
  static float nSigma(const TrackType& track, const o2::pid::tof::TOFResoParamsV3& params = parameters)
  {
    return nSigma<id>(track.tofSignal(), track.tofExpMom(), track.length(), track.p(), track.eta(), track.tofEvTime(), track.tofEvTimeErr(), params);
  }

  static bool isInit() { return mIsInit; } //! Get the initialization flag

  // Getters for the configurable options
  bool cfgAutoSetProcessFunctions() const { return mAutoSetProcessFunctions; }
  o2::common::core::CollisionSystemType::collType cfgCollisionType() const { return mCollisionSystem; }

 private:
  void inheritFromBaseTask(o2::framework::InitContext& initContext, const std::string task = "tof-signal");

  static bool mIsInit;       //! Flag to check if the parameters are initialized
  static int mLastRunNumber; //! Last run number for which the calibration was loaded

  o2::ccdb::BasicCCDBManager* mCcdb = nullptr; // Pointer to the CCDB manager

  // Configurable options
  std::string mUrl = "undefined";
  std::string mPathGrpLhcIf = "undefined";
  int64_t mTimestamp = -1;
  std::string mTimeShiftCCDBPathPos = "undefined";
  std::string mTimeShiftCCDBPathNeg = "undefined";
  std::string mTimeShiftCCDBPathPosMC = "undefined";
  std::string mTimeShiftCCDBPathNegMC = "undefined";
  std::string mParamFileName = "undefined";
  std::string mParametrizationPath = "undefined";
  std::string mReconstructionPass = "undefined";
  std::string mReconstructionPassDefault = "undefined";
  bool mFatalOnPassNotAvailable = false;
  bool mEnableTimeDependentResponse = false;
  o2::common::core::CollisionSystemType::collType mCollisionSystem = o2::common::core::CollisionSystemType::kCollSysUndef;
  bool mAutoSetProcessFunctions = false;

  template <typename VType>
  void getCfg(o2::framework::InitContext& initContext, const std::string name, VType& v, const std::string task)
  {
    if (!getTaskOptionValue(initContext, task, name, v, false)) {
      LOG(fatal) << "Could not get " << name << " from " << task << " task";
    }
  }
};

struct TOFResponse : o2::framework::LoadableServicePlugin<TOFResponseImpl> {
  TOFResponse() : LoadableServicePlugin{"O2PhysicsAnalysisCore:TOFSupport"}
  {
  }
};

} // namespace o2::pid::tof

#endif // COMMON_CORE_PID_PIDTOFPARAMSERVICE_H_
