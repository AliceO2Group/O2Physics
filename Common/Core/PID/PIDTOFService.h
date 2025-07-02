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
/// \file   PIDTOFService.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  30/06/2025
/// \brief  Implementation of the TOF PID service for the detector response
///

#ifndef COMMON_CORE_PID_PIDTOFSERVICE_H_
#define COMMON_CORE_PID_PIDTOFSERVICE_H_

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

#endif // COMMON_CORE_PID_PIDTOFSERVICE_H_
