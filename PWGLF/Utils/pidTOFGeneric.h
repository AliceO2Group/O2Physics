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
/// \file pidTOFGeneric.h
/// \brief Utilities to recalculate secondary tracks TOF PID
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
///

#ifndef PWGLF_UTILS_PIDTOFGENERIC_H_
#define PWGLF_UTILS_PIDTOFGENERIC_H_
#include "CollisionTypeHelper.h"
#include "MetadataHelper.h"
#include "TableHelper.h"

#include "Common/Core/PID/PIDTOF.h"

#include "CommonDataFormat/InteractionRecord.h"

#include <map>
#include <string>

namespace o2::aod
{

namespace pidtofgeneric
{

// Configuration common to all tasks, copied from pidTOFMerge.cxx but add metadataInfo as a member variable
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
  void getCfg(o2::framework::InitContext& initContext, const std::string name, VType& v, const std::string task)
  {
    if (!getTaskOptionValue(initContext, task, name, v, false)) {
      LOG(fatal) << "Could not get " << name << " from " << task << " task";
    }
  }

  void inheritFromBaseTask(o2::framework::InitContext& initContext, const std::string task = "tof-signal")
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
        mReconstructionPass = metadataInfo.get("RecoPassName");
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
            if (metadataInfo.isRun3()) {
              mRespParamsV3.setResolutionParametrization(paramCollection.getPars(mReconstructionPassDefault));
            } else {
              mRespParamsV3.setResolutionParametrizationRun2(paramCollection.getPars(mReconstructionPassDefault));
            }
            mRespParamsV3.setMomentumChargeShiftParameters(paramCollection.getPars(mReconstructionPassDefault));
          }
        }
      } else { // Pass is available, load non standard parameters
        if (metadataInfo.isRun3()) {
          mRespParamsV3.setResolutionParametrization(paramCollection.getPars(mReconstructionPass));
        } else {
          mRespParamsV3.setResolutionParametrizationRun2(paramCollection.getPars(mReconstructionPass));
        }
        mRespParamsV3.setMomentumChargeShiftParameters(paramCollection.getPars(mReconstructionPass));
      }
    } else if (!mEnableTimeDependentResponse) { // Loading it from CCDB
      LOG(info) << "Loading initial exp. sigma parametrization from CCDB, using path: " << mParametrizationPath << " for timestamp " << mTimestamp;
      o2::tof::ParameterCollection* paramCollection = ccdb->template getSpecific<o2::tof::ParameterCollection>(mParametrizationPath, mTimestamp);
      if (!paramCollection->retrieveParameters(mRespParamsV3, mReconstructionPass)) { // Attempt at loading the parameters with the pass defined
        if (mFatalOnPassNotAvailable) {
          LOG(fatal) << "Pass '" << mReconstructionPass << "' not available in the retrieved CCDB object";
        } else {
          LOG(warning) << "Pass '" << mReconstructionPass << "' not available in the retrieved CCDB object, fetching '" << mReconstructionPassDefault << "'";
          if (!paramCollection->retrieveParameters(mRespParamsV3, mReconstructionPassDefault)) {
            paramCollection->print();
            LOG(fatal) << "Cannot get default pass for calibration " << mReconstructionPassDefault;
          } else {
            if (metadataInfo.isRun3()) {
              mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPassDefault));
            } else {
              mRespParamsV3.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPassDefault));
            }
            mRespParamsV3.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPassDefault));
          }
        }
      } else { // Pass is available, load non standard parameters
        if (metadataInfo.isRun3()) {
          mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPass));
        } else {
          mRespParamsV3.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPass));
        }
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
        ccdb->setFatalWhenNull(false);
        mRespParamsV3.setTimeShiftParameters(ccdb->template getSpecific<TGraph>(nameShift, mTimestamp, metadata), isPositive);
        ccdb->setFatalWhenNull(true);
      }
      LOG(info) << " test getTimeShift at 0 " << (isPositive ? "pos" : "neg") << ": "
                << mRespParamsV3.getTimeShift(0, isPositive);
    };

    const std::string nameShiftPos = metadataInfo.isMC() ? mTimeShiftCCDBPathPosMC : mTimeShiftCCDBPathPos;
    updateTimeShift(nameShiftPos, true);
    const std::string nameShiftNeg = metadataInfo.isMC() ? mTimeShiftCCDBPathNegMC : mTimeShiftCCDBPathNeg;
    updateTimeShift(nameShiftNeg, false);

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
    LOG(info) << "Updating the parametrization from last run " << mLastRunNumber << " to " << bc.runNumber() << " and timestamp from " << mTimestamp << " " << bc.timestamp();
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
            if (metadataInfo.isRun3()) {
              mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPassDefault));
            } else {
              mRespParamsV3.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPassDefault));
            }
            mRespParamsV3.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPassDefault));
          }
        }
      } else { // Found the non default case
        if (metadataInfo.isRun3()) {
          mRespParamsV3.setResolutionParametrization(paramCollection->getPars(mReconstructionPass));
        } else {
          mRespParamsV3.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPass));
        }
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
      ccdb->setFatalWhenNull(false);
      mRespParamsV3.setTimeShiftParameters(ccdb->template getSpecific<TGraph>(nameShift, mTimestamp, metadata), isPositive);
      ccdb->setFatalWhenNull(true);
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

  o2::common::core::MetadataHelper metadataInfo; // additional member variable to store metadata information compared to pidTOFMerge.cxx

 private:
  int mLastRunNumber = -1; // Last run number for which the calibration was loaded
  int mInitMode = 0;       // 0: no init, 1: init, 2: inherit

  // Configurable options
  std::string mUrl;
  std::string mPathGrpLhcIf;
  int64_t mTimestamp;
  std::string mTimeShiftCCDBPathPos;
  std::string mTimeShiftCCDBPathNeg;
  std::string mTimeShiftCCDBPathPosMC;
  std::string mTimeShiftCCDBPathNegMC;
  std::string mParamFileName;
  std::string mParametrizationPath;
  std::string mReconstructionPass;
  std::string mReconstructionPassDefault;
  bool mFatalOnPassNotAvailable;
  bool mEnableTimeDependentResponse;
  int mCollisionSystem;
  bool mAutoSetProcessFunctions;
};

static constexpr float kCSPEED = TMath::C() * 1.0e2f * 1.0e-12f; // c in cm/ps

template <typename TTrack>
class TofPidNewCollision
{
 public:
  TofPidNewCollision() = default;
  ~TofPidNewCollision() = default;

  o2::pid::tof::TOFResoParamsV3 mRespParamsV3;
  o2::track::PID::ID pidType;

  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<TTrack, pid>;
  static constexpr auto responseEl = ResponseImplementation<o2::track::PID::Electron>();
  static constexpr auto responseMu = ResponseImplementation<o2::track::PID::Muon>();
  static constexpr auto responsePi = ResponseImplementation<o2::track::PID::Pion>();
  static constexpr auto responseKa = ResponseImplementation<o2::track::PID::Kaon>();
  static constexpr auto responsePr = ResponseImplementation<o2::track::PID::Proton>();
  static constexpr auto responseDe = ResponseImplementation<o2::track::PID::Deuteron>();
  static constexpr auto responseTr = ResponseImplementation<o2::track::PID::Triton>();
  static constexpr auto responseHe = ResponseImplementation<o2::track::PID::Helium3>();
  static constexpr auto responseAl = ResponseImplementation<o2::track::PID::Alpha>();

  void SetParams(o2::pid::tof::TOFResoParamsV3 const& para)
  {
    mRespParamsV3.setParameters(para);
  }

  void SetPidType(o2::track::PID::ID pidId)
  {
    pidType = pidId;
  }

  template <typename ParamType, typename TCollision>
  float GetTOFNSigma(o2::track::PID::ID pidId, const ParamType& parameters, TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D = true);

  template <typename ParamType, typename TCollision>
  float GetTOFNSigma(const ParamType& parameters, TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D = true);

  template <typename ParamType>
  float GetTOFNSigma(const ParamType& parameters, TTrack const& track);
  template <typename ParamType>
  float GetTOFNSigma(o2::track::PID::ID pidId, const ParamType& parameters, TTrack const& track);

  template <typename ParamType>
  float CalculateTOFNSigma(o2::track::PID::ID pidId, const ParamType& parameters, TTrack const& track, double tofsignal, double evTime, double evTimeErr)
  {

    float expSigma, tofNsigma = -999;

    switch (pidId) {
      case 0:
        expSigma = responseEl.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseEl.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 1:
        expSigma = responseMu.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseMu.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 2:
        expSigma = responsePi.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responsePi.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 3:
        expSigma = responseKa.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseKa.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 4:
        expSigma = responsePr.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responsePr.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 5:
        expSigma = responseDe.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseDe.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 6:
        expSigma = responseTr.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseTr.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 7:
        expSigma = responseHe.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseHe.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      case 8:
        expSigma = responseAl.GetExpectedSigma(parameters, track, tofsignal, evTimeErr);
        tofNsigma = (tofsignal - evTime - responseAl.GetCorrectedExpectedSignal(parameters, track)) / expSigma;
        break;
      default:
        LOG(fatal) << "Wrong particle ID in TofPidSecondary class";
        return -999;
    }

    return tofNsigma;
  }
};

template <typename TTrack>
template <typename ParamType, typename TCollision>
float TofPidNewCollision<TTrack>::GetTOFNSigma(o2::track::PID::ID pidId, const ParamType& parameters, TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D)
{

  if (!track.has_collision() || !track.hasTOF()) {
    return -999;
  }

  float mMassHyp = o2::track::pid_constants::sMasses2Z[track.pidForTracking()];
  float expTime = track.length() * std::sqrt((mMassHyp * mMassHyp) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v

  float evTime = correctedcol.evTime();
  float evTimeErr = correctedcol.evTimeErr();
  float tofsignal = track.trackTime() * 1000 + expTime; // in ps

  if (originalcol.globalIndex() == correctedcol.globalIndex()) {
    evTime = track.evTimeForTrack();
    evTimeErr = track.evTimeErrForTrack();
  } else {
    if (EnableBCAO2D) {
      auto originalbc = originalcol.template bc_as<o2::aod::BCsWithTimestamps>();
      auto correctedbc = correctedcol.template bc_as<o2::aod::BCsWithTimestamps>();
      o2::InteractionRecord originalIR = o2::InteractionRecord::long2IR(originalbc.globalBC());
      o2::InteractionRecord correctedIR = o2::InteractionRecord::long2IR(correctedbc.globalBC());
      tofsignal += originalIR.differenceInBCNS(correctedIR) * 1000;
    } else {
      auto originalbc = originalcol.template foundBC_as<o2::aod::BCsWithTimestamps>();
      auto correctedbc = correctedcol.template foundBC_as<o2::aod::BCsWithTimestamps>();
      o2::InteractionRecord originalIR = o2::InteractionRecord::long2IR(originalbc.globalBC());
      o2::InteractionRecord correctedIR = o2::InteractionRecord::long2IR(correctedbc.globalBC());
      tofsignal += originalIR.differenceInBCNS(correctedIR) * 1000;
    }
  }

  float tofNsigma = CalculateTOFNSigma(pidId, parameters, track, tofsignal, evTime, evTimeErr);
  return tofNsigma;
}

template <typename TTrack>
template <typename ParamType, typename TCollision>
float TofPidNewCollision<TTrack>::GetTOFNSigma(const ParamType& parameters, TTrack const& track, TCollision const& originalcol, TCollision const& correctedcol, bool EnableBCAO2D)
{
  return GetTOFNSigma(pidType, parameters, track, originalcol, correctedcol, EnableBCAO2D);
}

template <typename TTrack>
template <typename ParamType>
float TofPidNewCollision<TTrack>::GetTOFNSigma(o2::track::PID::ID pidId, const ParamType& parameters, TTrack const& track)
{

  if (!track.has_collision() || !track.hasTOF()) {
    return -999;
  }

  float mMassHyp = o2::track::pid_constants::sMasses2Z[track.pidForTracking()];
  float expTime = track.length() * std::sqrt((mMassHyp * mMassHyp) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v

  float evTime = track.evTimeForTrack();
  float evTimeErr = track.evTimeErrForTrack();
  float tofsignal = track.trackTime() * 1000 + expTime; // in ps

  float tofNsigma = CalculateTOFNSigma(pidId, parameters, track, tofsignal, evTime, evTimeErr);
  return tofNsigma;
}

template <typename TTrack>
template <typename ParamType>
float TofPidNewCollision<TTrack>::GetTOFNSigma(const ParamType& parameters, TTrack const& track)
{
  return GetTOFNSigma(pidType, parameters, track);
}

} // namespace pidtofgeneric
} // namespace o2::aod

#endif // PWGLF_UTILS_PIDTOFGENERIC_H_
